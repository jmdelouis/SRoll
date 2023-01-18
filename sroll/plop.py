import numpy as np
import sys
#import matplotlib.pyplot as plt
#import healpy as hp
import time
import tensorflow.compat.v1 as tf
import os
#test
from collections import defaultdict
import psutil as ps
from netCDF4 import Dataset



#Temporary global parameters
doloss2=False
nfile = 1
iocef = 0 # POIDS DES PWSTs
nout = 32 #nside map de sortie

tf.disable_v2_behavior()

if os.getenv('HOST')=='br146-050':
  lrank=int(os.getenv('OMPI_COMM_WORLD_NODE_RANK'))
else:
  lrank=0

os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"   # see issue #152
#TO BE DEBUGGED, NUMBER OF GPUS FORCED TO 3 SHOULD BE CLEANED
os.environ["CUDA_VISIBLE_DEVICES"]="%d"%(lrank%3)
print('RANK %d uses the GPU %d'%(lrank,lrank%3))
tf_device='/gpu:0'

#For run sur Garoupe
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # Only errors displayed

tf.logging.set_verbosity(tf.logging.ERROR)  # Double check: only errors displayed

#plt.switch_backend('agg') # Backend for plots without displaying them (OCCIGEN)

if lrank==0:
  print("Num GPUs Available: ", len(tf.config.experimental.list_physical_devices('GPU')))
tf.debugging.set_log_device_placement(False)


gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
  try:
    # Currently, memory growth needs to be the same across GPUs
    for gpu in gpus:
      tf.config.experimental.set_memory_growth(gpu, True)
    logical_gpus = tf.config.experimental.list_logical_devices('GPU')
    if lrank==0:
      print(len(gpus), "Physical GPUs,", len(logical_gpus), "Logical GPUs")
  except RuntimeError as e:
    # Memory growth must be set before GPUs have been initialized
    print(e)



# ----------------------------------------------------------------------------------------
def plt_vec(data):
  plt.plot(data)
  plt.show()
  del(data)
# ---------------------------------------------------------------------------------------- 
def plt_histo(data,hdata):
  idx=((data*256)/(2*np.pi)).astype('int')
  idx[idx==256]=255

  hh=np.bincount(idx,weights=hdata,minlength=256)

  if hh.std()>200:
    print('VAR %.4g %.4g'%(hh.std(),hh.mean()))
    plt.plot(data[hdata>0])
    plt.show()
    plt.plot(hh)
    plt.show()
  del(data)
  del(hdata)
# ----------------------------------------------------------------------------------------
def lrelu(x, alpha=0.3):
  #return tf.nn.relu(x)
  return tf.maximum(x, tf.multiply(x, alpha))
# ---------------------------------------------------------------------------------------- 
def wdecod(parameter,netinfo,rank): 
  """convert into image/bias"""
 
  #print(parameter.get_shape().as_list())
  if rank==0:
    sys.stderr.write('NPAR PER IMAGE %d\n'%(int((parameter.get_shape().as_list())[0])))
  if netinfo['DOFULLYCONNECTED']==False:
    hidden = parameter
  else:
    hidden = tf.matmul(parameter, netinfo['fc2_weights']) + netinfo['fc2_biases']
    
  if netinfo['DORELU']==True:
      hidden = lrelu(hidden)
      #hidden = lrelu(hidden)

  if rank==0:
    sys.stderr.write('NPAR SIZE %d %d\n'%(int((parameter.get_shape().as_list())[0]),int((parameter.get_shape().as_list())[1])))
  NTIME=int((hidden.get_shape().as_list())[0])

  relu = tf.reshape(hidden,[NTIME, netinfo['XIMAGE_SIZE'] // (netinfo['XSCALE']**netinfo['NDCONV']), netinfo['YIMAGE_SIZE']// (netinfo['YSCALE']**(netinfo['NDCONV'])), netinfo['NUM_DCONV_CHAN'][0]])


  ## deconvolve : First layer
  # DeConvolution kernel



  for i in range(netinfo['NDCONV']):
   
    conv = tf.nn.conv2d_transpose(relu,netinfo['dconv_weights'][i],strides=[1, netinfo['XSCALE'],netinfo['YSCALE'] , 1],padding='SAME', 
                                  output_shape=[NTIME,netinfo['XIMAGE_SIZE'] // (netinfo['XSCALE']**(netinfo['NDCONV']-1-i)),
                                  netinfo['YIMAGE_SIZE']// (netinfo['YSCALE']**(netinfo['NDCONV']-1-i)),netinfo['NUM_DCONV_CHAN'][i+1]],name='dconv_%d'%(i))
    #Print for debug    
    #print(conv.get_shape().as_list())
    #print(netinfo['dconv_biases'][i].get_shape().as_list())
    #print(netinfo['dconv_weights'][i].get_shape().as_list())
    
    tmp=tf.nn.bias_add(conv, netinfo['dconv_biases'][i])
    if rank==0:
      alist=tmp.get_shape().as_list()
      s="["
      for iv in alist:
        s=s+"%d,"%(iv)
      s=s+"],["
      alist=netinfo['dconv_weights'][i].get_shape().as_list()
      for iv in alist:
        s=s+"%d,"%(iv)
      s=s+"]\n"
      sys.stderr.write(s)

    if netinfo['DORELU']==True:
      if i!=netinfo['NDCONV']-1:
        relu = lrelu(tmp)
      else:
        #relu = lrelu(tmp,alpha=0)
        relu=tmp
    else:
      relu=tmp

  relu=relu/100

  return relu
# ---------------------------------------------------------------------------------------- 
def cnn_decod(parameter,map_weights,map_biases):

    param = parameter
    if DORELU==True:
        param = lrelu(param)

    relu_shape = param.get_shape().as_list()
    print(relu_shape,param.get_shape().as_list())
    
    relu=param
    # DeConvolution kernel
    for i in range(mapdepth):
        relu_shape = relu.get_shape().as_list()
        outchan=mapchan
        if i==mapdepth-1:
            outchan=2
        conv = tf.nn.conv2d_transpose(relu,map_weights[i],strides=[1, 4, 1, 1],padding='SAME',
                                      output_shape=[relu_shape[0],4*relu_shape[1],1,outchan])

        relu=tf.nn.bias_add(conv, map_biases[i])
        if DORELU==True and i!=mapdepth-1:
            relu = lrelu(relu)
            print('NL TRANSFORM')
        relu_shape = relu.get_shape().as_list()

        print(relu_shape)

    relu=tf.reshape(relu,[relu_shape[0]*relu_shape[1],relu_shape[3]])/100
    return relu
# ----------------------------------------------------------------------------------------
def model(parameter,netinfo,rank):
  modpar=wdecod(parameter,netinfo,rank)
  return modpar
# ----------------------------------------------------------------------------------------
def hpwst(image1,image2,mask,first=True,doheavy=False):
    """Calcul des Pwst avec option heavy"""

    #Init
    im_shape = image1.get_shape().as_list()
    nout=int(im_shape[1])
    nstep=int(np.log(np.sqrt(nout/12))/np.log(2))
    lim1=image1
    lim2=image2
    tshape=mask.get_shape().as_list()

    
    if doheavy==True:
        npar=8
    else:
        npar=1
    if first==True:
        vmask=tf.tile(mask,tf.constant([1,1,1,4], tf.int32))
    else:
        vmask=tf.tile(mask,tf.constant([1,1,1,npar*4], tf.int32))
    lim1=image1
    lim2=image2
    
    n0=int(np.sqrt(nout/12))
    for iscale in range(nstep):

        if first==True:
            vnorm=4.0/(tf.math.reduce_sum(vmask))
        else:
            vnorm=npar*16.0/(tf.math.reduce_sum(vmask))

        im_shape = lim1.get_shape().as_list()
        #print(np.abs(np.mean(np.dot(map.reshape(6,8),ww).reshape(6,4,2),1)))
        if first==True:
            lconv1 = lim1 - tf.reshape(tf.reduce_mean(tf.reshape(tf.gather(lim1, widx[n0],axis=1),[1,12*n0*n0,4,1,2]),2),[1,12*n0*n0,1,2])
            lconv2 = lim2 - tf.reshape(tf.reduce_mean(tf.reshape(tf.gather(lim2, widx[n0],axis=1),[1,12*n0*n0,4,1,2]),2),[1,12*n0*n0,1,2])
            conv1= tf.reshape(tf.tile(lconv1,[1,1,4,1]),[1,12*n0*n0,1,8]) + tf.reshape(tf.gather(lconv1, widx2[n0],axis=1),[1,12*n0*n0,1,8])
            conv2= tf.reshape(tf.tile(lconv2,[1,1,4,1]),[1,12*n0*n0,1,8]) + tf.reshape(tf.gather(lconv2, widx2[n0],axis=1),[1,12*n0*n0,1,8])
            tconv=conv1*conv2
            tmp=tf.sign(tconv)*tf.sqrt(tf.sign(tconv)*tconv)
        else:
            conv1 = lim1 - tf.reduce_mean(tf.reshape(tf.gather(lim1, widx[n0],axis=1),[1,12*n0*n0,4,1,8]),2)
            conv1 = tf.reshape(tf.tile(conv1,[1,1,4,1]),[1,12*n0*n0,1,32]) + tf.reshape(tf.gather(conv1, widx2[n0],axis=1),[1,12*n0*n0,1,32])
            if doheavy==True:
                conv1=tf.reshape(conv1,[1,12*n0*n0,1,4,8])
                tconv=tf.reshape(tf.reshape(tf.repeat(conv1,8),[1,12*n0*n0,1,4,8,8])*tf.reshape(tf.tile(conv1,[1,1,1,1,8]),[1,12*n0*n0,1,4,8,8]),[1,12*n0*n0,1,256])
                tmp=tf.sign(tconv)*tf.sqrt(tf.sign(tconv)*tconv)
            else:
                tmp=tf.abs(conv1) #tf.sqrt(tf.sign(tconv)*tconv)

        #tmp=tf.nn.avg_pool(tmp, ksize=[1, 2, 2, 1], strides=[1, 2, 2, 1],padding='SAME')
        val=vnorm*tf.math.reduce_sum(vmask*tmp,1)
        tshape=tmp.get_shape().as_list()
        if iscale==0:
            s1=tf.reshape(val,[1,4,tshape[3]//4])
        else:
            s1=tf.concat([s1,tf.reshape(val,[1,4,tshape[3]//4])],0)

        if first==True and iscale<nstep-1:
            val2=(hpwst(tf.nn.relu(tmp),tf.nn.relu(tmp),vmask,first=False,doheavy=doheavy)-hpwst(tf.nn.relu(-tmp),tf.nn.relu(-tmp),vmask,first=False,doheavy=doheavy))
            if iscale==0:
                i1=tf.constant(np.repeat(np.zeros(nstep,dtype='int'),npar*32).reshape(nstep,npar*4,4,2))
                i2=tf.constant(np.repeat(np.arange(nstep,dtype='int'),npar*32).reshape(nstep,npar*4,4,2))
                s2=tf.reshape(val2,[nstep,npar*4,4,2])
            else:
                s2=tf.concat([s2,tf.reshape(val2,[nstep-iscale,npar*4,4,2])],0)
                i1=tf.concat([i1,tf.constant(np.repeat(iscale+np.zeros(nstep-iscale,dtype='int'),npar*32).reshape(nstep-iscale,npar*4,4,2))],0)
                i2=tf.concat([i2,tf.constant(np.repeat(iscale+np.arange(nstep-iscale,dtype='int'),npar*32).reshape(nstep-iscale,npar*4,4,2))],0)

        lim1=0.5*tf.nn.avg_pool(lim1, ksize=[1, 4, 1, 1], strides=[1, 4, 1, 1],padding='SAME')
        lim2=0.5*tf.nn.avg_pool(lim2, ksize=[1, 4, 1, 1], strides=[1, 4, 1, 1],padding='SAME')
        vmask=tf.nn.avg_pool(vmask, ksize=[1, 4, 1, 1], strides=[1, 4, 1, 1],padding='SAME')
        n0=n0//2

    if first==True:
        return(s1,s2,i1,i2)
    else:
        return(s1)
#----------------------------------------------------------------------------------------
def init_calcul_pwst(TFLEARN,nout):
    """ Inititalisation calcul des PWST 
        Params : 
            - nout = nside map de sortie
        Return :
            widx = tab[tf.constant]
            widx2 = tab[tf.constant]
    """
    #Call in init_nework
    #Temporay defined here : 
    #PARAM    
    path = TFLEARN+'/WIDXR_%d.npy'    
    nstep=int(np.log(nout)/np.log(2))+1
    widx={}
    widx2={}
    for i in range(nstep):
        lout=nout//(2**i)
        tmp=((np.load(path%(lout)).flatten()).astype('int32')).reshape(12*lout**2,8)
        widx[lout]=tf.constant(tmp[:,0:4])
        widx2[lout]=tf.constant(tmp[:,[0,4,1,7]])

        #print('WAVELET INDEX',lout)

    return widx,widx2
# ----------------------------------------------------------------------------------------
def down(im,nside):
    nin=int(np.sqrt(im.shape[0]/12))
    if nin==nside:
        return(im)
    nn=(nin//nside)**2
    return(np.mean(im.reshape(12*nin**2//nn,nn),1))
# ----------------------------------------------------------------------------------------
def alloc_table_float32(ndata):
  mydata = np.zeros([ndata],dtype='float32')
  return mydata
# ---------------------------------------------------------------------------------------- 
def alloc_table_int32(ndata):
  mydata = np.zeros([ndata],dtype='int32')
  return mydata
# ---------------------------------------------------------------------------------------- 
def free_table(data):
  del(data)
# ---------------------------------------------------------------------------------------- 
def init_shape(nbolo,xsize,ysize,rank,
               ncomp,ndconv,kernelsz,scale,npar_fc,
               loss_choice,flag_normalize_data_cnn_loss,huber_loss_delta,
               num_epochs,eval_freq,learning_rate,decay_rate,batchsz,
               flag_do_transfer_learning,flag_do_fully_connected,flag_do_relu,
               flag_test_convergence,flag_no_net,
               flag_save_cnn_weights, path_cnn_weights,
               seed):
  """Init shape for CNN  """
  

  mynetwork = {}
  mynetwork['FromFile'] = False

  # Network dimensions
  mynetwork['nbolo'] = nbolo
  mynetwork['XIMAGE_SIZE'] = xsize
  mynetwork['YIMAGE_SIZE'] = ysize
  mynetwork['rank'] = rank
  TESTCONVERGENCE = False
  NH = 4
  nonet = False
  # Network architecture
  mynetwork['NCOMP'] = ncomp
  mynetwork['NDCONV'] = ndconv
  mynetwork['KERNELSZ'] = kernelsz
  mynetwork['SCALE'] = scale
  mynetwork['DOFULLYCONNECTED'] = flag_do_fully_connected
  mynetwork['DORELU'] = flag_do_relu
  
  XNHIDDEN = mynetwork['XIMAGE_SIZE'] // (mynetwork['SCALE']**mynetwork['NDCONV'])
  YNHIDDEN = mynetwork['YIMAGE_SIZE'] 
  #NHIDDEN = mynetwork['NCOMP']
  DEEPNESS = mynetwork['NCOMP']
     
  if mynetwork['DOFULLYCONNECTED']:
      NHIDDEN = NH
  else:
      NHIDDEN = mynetwork['NCOMP']

  if mynetwork['DOFULLYCONNECTED']:
      NHIDDEN = npar_fc
  else:
      NHIDDEN = mynetwork['NCOMP']
  
  NUM_DCONV_CHAN = {}



  #INIT_learning  -------------------------------------
  # CNN loss function parameters
  mynetwork['CNN_LOSS'] = loss_choice   
  mynetwork['CNN_NORMALIZE_DATA_IN_LOSS'] = flag_normalize_data_cnn_loss
  mynetwork['CNN_HUBER_LOSS_DELTA'] = huber_loss_delta
               

  # Training parameters
  mynetwork['NUM_EPOCHS'] = num_epochs
  EVAL_FREQUENCY = eval_freq
  mynetwork['LEARNING_RATE']= learning_rate
  mynetwork['DECAY_RATE']=decay_rate
  mynetwork['BATCHSZ'] = batchsz
  # --------------------------------------------------------



  # Debug options
  TESTCONVERGENCE = flag_test_convergence
  nonet = flag_no_net
  
  # Transfer learning
  mynetwork['DOTRANSFERLEARNING'] = flag_do_transfer_learning

  # Save weights parametres
  mynetwork['SAVE_CNN_WEIGHTS']=flag_save_cnn_weights
  mynetwork['CNN_WEIGHTS'] = path_cnn_weights

  # Random number generator seed
  SEED = seed
  #=================================================================================================
  #   NEXT LINES CAN BE USED TO CHECK THE CONVERGENCE
  #=================================================================================================
  if TESTCONVERGENCE==True:
    model2=np.zeros([8,256,256])
    xx=np.arange(256,dtype='float')
    model2[0,:,:]=np.exp(-1E-3*((np.repeat(xx,256)-128)**2+(np.tile(xx,256)-140)**2)).reshape(256,256)
    model2-=np.mean(model2)
    model2=model2.flatten()
    np.save('model2.npy',model2)
    signal=(0.1*np.random.randn(signal.shape[0])+model2[idx]).astype('float32')

    x,y,z=hp.pix2vec(2048,np.arange(12*2048*2048))

    signal+=x[hidx]+TCO1*y[hidx]+TSI1*z[hidx]


  if mynetwork['DOFULLYCONNECTED']==False:

    NTOTHIDDEN=NHIDDEN*(XNHIDDEN*YNHIDDEN)
    
    for i in range(mynetwork['NDCONV']):
      NUM_DCONV_CHAN[i] = int(DEEPNESS)

    NUM_DCONV_CHAN[0] = NHIDDEN

  else:
    NTOTHIDDEN=NHIDDEN
    for i in range(mynetwork['NDCONV']):
      NUM_DCONV_CHAN[i]=DEEPNESS

  NUM_DCONV_CHAN[mynetwork['NDCONV']]=1

  mynetwork['NUM_DCONV_CHAN']=NUM_DCONV_CHAN
  mynetwork['NTOTHIDDEN']=NTOTHIDDEN
  mynetwork['RAPDATA']=1.0

  if mynetwork['rank']==0:
    sys.stderr.write('init shape done\n')

  return(mynetwork)
# ----------------------------------------------------------------------------------------
def init_shape_from_files(TFLEARN,bolotab,rank,CNN_LEARN_MODEL,tmpname):
    """ 
    Init shape of CNNs FSL and MAP form files  

    => Parameters : 			  
            TFLEARN :  tab of NET_INFO files -- CNN FSL/MAP
            bolotab : taf of string -- names of bolometers
            rank : ?
            tmpid : ?
            tmpname : Name used for the output (TEMPORARY DEBUG)

    => Return :
            my_network : dict 2Keys - mynetwork[][] 
    """   
    
    NFC        = {}
    LEARNPARAM = {}
    LEARNFL    = {}
    LEARNCONV  = {}
    
    NFC['MAP']        = ((CNN_LEARN_MODEL%10)//8)%2
    LEARNPARAM['MAP'] = ((CNN_LEARN_MODEL%10)//4)%2
    LEARNFL['MAP']    = ((CNN_LEARN_MODEL%10)//2)%2
    LEARNCONV['MAP']  = (CNN_LEARN_MODEL%10)%2
    
    NFC['FSL']        = ((CNN_LEARN_MODEL//10)//8)%2
    LEARNPARAM['FSL'] = ((CNN_LEARN_MODEL//10)//4)%2
    LEARNFL['FSL']    = ((CNN_LEARN_MODEL//10)//2)%2
    LEARNCONV['FSL']  = (CNN_LEARN_MODEL//10)%2

    ##Initialisation 
    #init network
    mynetwork = defaultdict(dict)
    #For debug temporary path defined here
    path_npy = TFLEARN

    if rank==0:
      sys.stderr.write('A MEM %.2fGB AVAIL %.2fGB\n'%(ps.virtual_memory().used/1024/1024/1024,ps.virtual_memory().available/1024/1024/1024))

    #init names for dict
    CNN_names = ["MAP","FSL"]

    #load netcdf 
    #netcdf = Dataset('/export/home/tfoulquier/python_scripts/CNN_netcdf.nc','r', format='NETCDF4')
    for elt in CNN_names:
      #load netcdf
      if(elt == 'MAP'):
        netcdf = Dataset('MAP_CNN_netcdf.nc','r', format='NETCDF4')
      else:
        netcdf = Dataset('FSL_CNN_netcdf.nc','r', format='NETCDF4')
      
      """
      #Load file
      #INFO=np.load(path_npy+'/CNN_%s_NET_INFO.npy'%(elt))
      INFO=np.load('/export/home/tfoulquier/CNN_%s_NET_INFO.npy'%(elt))
           
      # Read net_info 
      nbolo=int(INFO[0])
      xsize=int(INFO[1])
      ysize=int(INFO[2])
      NDCONV=int(INFO[3])
      DOFULLYCONNECTED=int(INFO[4])==1
     
      """

      nbolo=netcdf[elt].dimensions['batchsz'].size
      xsize=netcdf[elt].dimensions['ximage_size'].size
      ysize=netcdf[elt].dimensions['yimage_size'].size
      NDCONV=netcdf[elt].dimensions['ndconv'].size
      DOFULLYCONNECTED=(netcdf[elt].dimensions['dofullyconnected'].size)==1      

      """
      print(" == DEBUG :")
      print("NETINFO",int(INFO[0]),int(INFO[1]),int(INFO[2]),int(INFO[3]),int(INFO[4]))
      if(netcdf[elt].dimensions['batchsz'].size == int(INFO[0])):
        print("PLOP ***************************")
        print(type(netcdf[elt].dimensions['batchsz'].size))

      
      print("NETCDF",
      int(netcdf[elt].dimensions['batchsz'].size),
      netcdf[elt].dimensions['ximage_size'].size,
      netcdf[elt].dimensions['yimage_size'].size,
      netcdf[elt].dimensions['ndconv'].size,
      netcdf[elt].dimensions['dofullyconnected'].size)
      
      #print("##### "+str(elt)+" --  NDCONV =",int(INFO[3]))
      """
      
      mynetwork[elt]['FromFile']=True
      
      mynetwork[elt]['NFC']        = NFC[elt]
      mynetwork[elt]['LEARNPARAM'] = LEARNPARAM[elt]
      mynetwork[elt]['LEARNFL']    = LEARNFL[elt]
      mynetwork[elt]['LEARNCONV']  = LEARNCONV[elt]
      
      if NFC[elt]==False:
        if LEARNFL[elt] == False:
          mynetwork[elt]['fc2_weights'] = tf.constant(netcdf[elt].variables['fc2_weights'][:])
          mynetwork[elt]['fc2_biases'] = tf.constant(netcdf[elt].variables['fc2_biases'][:])
        else:
          mynetwork[elt]['fc2_weights'] = tf.Variable(netcdf[elt].variables['fc2_weights'][:])
          mynetwork[elt]['fc2_biases'] = tf.Variable(netcdf[elt].variables['fc2_biases'][:])
      
      mynetwork[elt]['dconv_weights']= {}
      mynetwork[elt]['dconv_biases']= {}

      print("DCONV layer",elt)
      for i in range(NDCONV):
        if LEARNCONV[elt] == False:
          mynetwork[elt]['dconv_weights'][i]=tf.constant(netcdf[elt].variables['dconv_weights_'+str(i)][:])
          mynetwork[elt]['dconv_biases'][i]=tf.constant(netcdf[elt].variables['dconv_biases_'+str(i)][:])
        else:
          mynetwork[elt]['dconv_weights'][i]=tf.Variable(netcdf[elt].variables['dconv_weights_'+str(i)][:])
          mynetwork[elt]['dconv_biases'][i]=tf.Variable(netcdf[elt].variables['dconv_biases_'+str(i)][:])         
      print("DCONV layer pass",elt)
      

      mynetwork[elt]['XIMAGE_SIZE']=xsize
      mynetwork[elt]['YIMAGE_SIZE']=ysize
      mynetwork[elt]['nbolo']=nbolo
      mynetwork[elt]['rank']=rank
      #mynetwork[elt]['NCOMP']=int(INFO[6])
      mynetwork[elt]['NCOMP']= netcdf[elt].dimensions['ncomp'].size
      #mynetwork[elt]['XSCALE']=int(INFO[7])
      mynetwork[elt]['XSCALE'] = netcdf[elt].dimensions['"scalex"'].size
      #mynetwork[elt]['YSCALE']=int(INFO[9])
      mynetwork[elt]['YSCALE'] = netcdf[elt].dimensions['scaley'].size
      #mynetwork[elt]['RAP']=(np.sqrt(1.0/INFO[11])).astype('float32')
      mynetwork[elt]['RAP']=(np.sqrt(1.0/(netcdf[elt].dimensions['corrnorm'].size).astype('float32')
      mynetwork[elt]['NDCONV']=NDCONV
      mynetwork[elt]['TMPNAME']=tmpname

      XNHIDDEN=mynetwork[elt]['XIMAGE_SIZE'] // (mynetwork[elt]['XSCALE']**mynetwork[elt]['NDCONV'])
      YNHIDDEN=mynetwork[elt]['YIMAGE_SIZE'] // (mynetwork[elt]['YSCALE']**mynetwork[elt]['NDCONV'])
      
      mynetwork[elt]['XNHIDDEN']=XNHIDDEN
      mynetwork[elt]['YNHIDDEN']=YNHIDDEN

      #NHIDDEN = int(INFO[5])
      NHIDDEN = netcdf[elt].dimensions['deepness'].size
      DEEPNESS = netcdf[elt].dimensions['ncomp'].size
      EVAL_FREQUENCY = 100
      
      #mynetwork[elt]['XKERNELSZ'] = int(INFO[8])
      #mynetwork[elt]['YKERNELSZ'] = int(INFO[10])
      mynetwork[elt]['XKERNELSZ'] = netcdf[elt].dimensions['xkernelsz'].size
      mynetwork[elt]['YKERNELSZ'] = netcdf[elt].dimensions['ykernelsz'].size

      mynetwork[elt]['NUM_EPOCHS']= 2000
      SEED = 1234
      
      mynetwork[elt]['DORELU'] = True
      mynetwork[elt]['LEARNING_RATE']= 0.03
      mynetwork[elt]['DECAY_RATE']=0.995
      mynetwork[elt]['DOFULLYCONNECTED'] = NFC[elt]==0
      mynetwork[elt]['DOTRANSFERLEARNING'] = True
          
      if elt == 'MAP':
        mynetwork[elt]['BATCHSZ'] = 1
      else:
        mynetwork[elt]['BATCHSZ'] = nbolo
      
      TESTCONVERGENCE = False

      NUM_DCONV_CHAN = {}
      
      if mynetwork[elt]['DOFULLYCONNECTED']==False:
        NTOTHIDDEN=NHIDDEN*(XNHIDDEN*YNHIDDEN)
        for i in range(mynetwork[elt]['NDCONV']):
          #NUM_DCONV_CHAN[i]=int(DEEPNESS*2**i)
          NUM_DCONV_CHAN[i] = int(DEEPNESS)
        NUM_DCONV_CHAN[0] = NHIDDEN

      else:
        NTOTHIDDEN=NHIDDEN
        for i in range(mynetwork[elt]['NDCONV']):
          #NUM_DCONV_CHAN[i]=DEEPNESS*2**i
          NUM_DCONV_CHAN[i]=DEEPNESS
        NUM_DCONV_CHAN[0] =NHIDDEN
      
      if elt == 'MAP':
        NUM_DCONV_CHAN[mynetwork[elt]['NDCONV']]=2 #DIFFERENT!
      else:
        NUM_DCONV_CHAN[mynetwork[elt]['NDCONV']]=1
      
      mynetwork[elt]['NUM_DCONV_CHAN']=NUM_DCONV_CHAN
      mynetwork[elt]['NTOTHIDDEN']=NTOTHIDDEN
      mynetwork[elt]['TFLEARN']=TFLEARN
      
    
    if mynetwork[elt]['rank']==0:
        sys.stderr.write('init shape done\n')
    
    
    return(mynetwork)

# ---------------------------------------------------------------------------------------- 

def alloc_param(mynetwork): 
  """ Ecrit en dur for debug -- a reecrire """
 
  param = dict()

  #for elt in mynetwork:
  param=np.zeros([mynetwork['MAP']['BATCHSZ']*mynetwork['MAP']['NTOTHIDDEN'] + mynetwork['FSL']['BATCHSZ']*mynetwork['FSL']['NTOTHIDDEN']],dtype='float32')
  #print("param = "+str(param[elt])+"  elt = "+str(elt)+"  lenght = "+str(len(param))+"\n")
    

  return(param)
# ---------------------------------------------------------------------------------------- 
def get_param(myrun):
	sess=myrun['sess']
	params=sess.run([myrun['params']])[0]
	return(params)
# ---------------------------------------------------------------------------------------- 
def alloc_convw(mynetwork):
  outobj={}
  for elt in mynetwork:
    for i in range(mynetwork[elt]['NDCONV']):
      outobj["%03d"%(i)] = np.zeros([mynetwork[elt]['XKERNELSZ'] *mynetwork[elt]['YKERNELSZ'] *mynetwork[elt]['NUM_DCONV_CHAN'][i+1]*mynetwork[elt]['NUM_DCONV_CHAN'][i]],dtype='float32')
  return(outobj)
# ---------------------------------------------------------------------------------------- 
def alloc_convb(mynetwork):
  outobj={}
  for elt in mynetwork:
    for i in range(mynetwork[elt]['NDCONV']):
      outobj["%03d"%(i)] = np.zeros([mynetwork[elt]['NUM_DCONV_CHAN'][i+1]],dtype='float32')
  return(outobj)
# ---------------------------------------------------------------------------------------- 
def get_convw(myrun):
  sess=myrun['sess']
  outobj={}
  for i in range(mynetwork['NDCONV']):
    outobj["%03d"%(i)] = sess.run([myrun['dconvw']])[0]
  return(outobj)
# ---------------------------------------------------------------------------------------- 
def get_convb(myrun):
  sess=myrun['sess']
  outobj={}
  for i in range(mynetwork['NDCONV']):
    outobj["%03d"%(i)] = sess.run([myrun['dconvb']])[0]
  return(outobj)
# ---------------------------------------------------------------------------------------- 
def alloc_flw(mynetwork):
  flw={}
  ib=0
  for elt in mynetwork:
    flw['%03d'%(ib)]=np.zeros([mynetwork[elt]['NTOTHIDDEN']* mynetwork[elt]['XIMAGE_SIZE'] // (mynetwork[elt]['XSCALE']**mynetwork[elt]['NDCONV']) * mynetwork[elt]['YIMAGE_SIZE']// (mynetwork[elt]['YSCALE']**mynetwork[elt]['NDCONV']) * mynetwork[elt]['NUM_DCONV_CHAN'][0]],dtype='float32')
  return(flw)
# ----------------------------------------------------------------------------------------
def alloc_flb(mynetwork):
  flb={}
  ib=0
  for elt in mynetwork:
    flb['%03d'%(ib)]=np.zeros([mynetwork[elt]['XIMAGE_SIZE'] // (mynetwork[elt]['XSCALE']**mynetwork[elt]['NDCONV']) * mynetwork[elt]['YIMAGE_SIZE']// (mynetwork[elt]['YSCALE']**mynetwork[elt]['NDCONV']) * mynetwork[elt]['NUM_DCONV_CHAN'][0]],dtype='float32')
  return(flb)
# ---------------------------------------------------------------------------------------- 
def get_flw(myrun):
  sess=myrun['sess']
  flw= sess.run([myrun['fc2_weights']])[0]
  return(flw)
# ---------------------------------------------------------------------------------------- 
def get_flb(myrun):
  sess=myrun['sess']
  flb= sess.run([myrun['fc2_biases']])[0]
  return(flb)
# ---------------------------------------------------------------------------------------- 
def init_network(mynetwork,signal,weights,TCO1,TSI1,hidx,idx,realpix,coef): #,in_param,in_convw,in_convb,in_flw,in_flb):
    nonet = False
    netinfo={}
    nmem=0.0
    print("INIT NETWORK -----")


    for i in mynetwork:
      netinfo[i]=mynetwork[i]
   
    if mynetwork["MAP"]['rank']==0:
      sys.stderr.write('B MEM %.2fGB AVAIL %.2fGB\n'%(ps.virtual_memory().used/1024/1024/1024,ps.virtual_memory().available/1024/1024/1024))

    NPIX=realpix.shape[0]

    coefregul=weights.sum()

    params ={}
    npar = 0
    for elt in netinfo:
      print("elt = ",elt)
      print("mynetwork[elt]['BATCHSZ'] ",mynetwork[elt]['BATCHSZ'])
      print("mynetwork[elt]['NTOTHIDDEN'] ",mynetwork[elt]['NTOTHIDDEN'])


        
      if mynetwork[elt]['NFC']==True:
        params[elt] = tf.Variable(np.zeros([mynetwork[elt]['BATCHSZ'],mynetwork[elt]['XIMAGE_SIZE'] // (mynetwork[elt]['XSCALE']**mynetwork[elt]['NDCONV'])*mynetwork[elt]['YIMAGE_SIZE']// (mynetwork[elt]['YSCALE']**(mynetwork[elt]['NDCONV']))*mynetwork[elt]['NUM_DCONV_CHAN'][0]],dtype='float32'),name='LEVEL0_%s'%(elt))

      else: 
        if mynetwork[elt]['LEARNPARAM']==True:
          if elt=='MAP':
            params[elt] = tf.Variable(np.load(mynetwork['MAP']['TFLEARN']+'/MAPPAR.npy').reshape(mynetwork[elt]['BATCHSZ'],mynetwork[elt]['NTOTHIDDEN']),name='PARAM_%s'%(elt))
          else:
            params[elt] = tf.Variable(np.load(mynetwork['MAP']['TFLEARN']+'/FSLPAR.npy').reshape(mynetwork[elt]['BATCHSZ'],mynetwork[elt]['NTOTHIDDEN']),name='PARAM_%s'%(elt))
            
        else:
          if elt=='MAP':
            plop = np.load(mynetwork['MAP']['TFLEARN']+'/MAPPAR.npy')
            
            params[elt] = tf.constant(np.load(mynetwork['MAP']['TFLEARN']+'/MAPPAR.npy').reshape(mynetwork[elt]['BATCHSZ'],mynetwork[elt]['NTOTHIDDEN']))
          else:
            params[elt] = tf.constant(np.load(mynetwork['MAP']['TFLEARN']+'/FSLPAR.npy').reshape(mynetwork[elt]['BATCHSZ'],mynetwork[elt]['NTOTHIDDEN']))
            
        npar+=mynetwork[elt]['BATCHSZ']*mynetwork[elt]['NTOTHIDDEN']
    
    
    if mynetwork["MAP"]['rank']==0:
      sys.stderr.write('C MEM %.2fGB AVAIL %.2fGB\n'%(ps.virtual_memory().used/1024/1024/1024,ps.virtual_memory().available/1024/1024/1024))
    
    ii_signal=np.bincount(hidx,weights=weights*signal,minlength=realpix.shape[0]).astype('float32')
    iiw_signal=np.bincount(hidx,weights=weights,minlength=realpix.shape[0]).astype('float32')
    iico1 = np.bincount(hidx,weights=weights*TCO1,minlength=realpix.shape[0]).astype('float32')
    iisi1 = np.bincount(hidx,weights=weights*TSI1,minlength=realpix.shape[0]).astype('float32')
    
    lidx=np.where(iiw_signal>0)[0]
    
    if mynetwork["MAP"]['rank']==0:
      sys.stderr.write('D MEM %.2fGB AVAIL %.2fGB\n'%(ps.virtual_memory().used/1024/1024/1024,ps.virtual_memory().available/1024/1024/1024))
        
    ii_signal[lidx]/=iiw_signal[lidx]
    iico1[lidx]/=iiw_signal[lidx]
    iisi1[lidx]/=iiw_signal[lidx]

    iiw_signal[lidx] = 1 /iiw_signal[lidx]
    iiw_signal[iiw_signal<=0]=0

    signal=signal-ii_signal[hidx]
    
    del(ii_signal)
    del(lidx)
    
    bidx=idx//(netinfo["FSL"]['XIMAGE_SIZE']*netinfo["FSL"]['YIMAGE_SIZE'])
    OUT_CO = TCO1
    OUT_SI = TSI1
    CO = (TCO1 - iico1[hidx])
    SI = (TSI1 - iisi1[hidx])
    
    del(bidx)
    del(TCO1)
    del(TSI1)
    del(iisi1)
    del(iico1)
    
    if mynetwork["MAP"]['rank']==0:
      sys.stderr.write('E MEM %.2fGB AVAIL %.2fGB\n'%(ps.virtual_memory().used/1024/1024/1024,ps.virtual_memory().available/1024/1024/1024))
    
    tf_CO  = tf.constant(CO)
    nmem+=np.array(CO.shape).sum()
    tf_SI  = tf.constant(SI)
    nmem+=np.array(SI.shape).sum()
    tf_OUT_CO = tf.constant(OUT_CO)
    nmem+=np.array(OUT_CO.shape).sum()
    tf_OUT_SI = tf.constant(OUT_SI)
    nmem+=np.array(OUT_SI.shape).sum()

    data   = tf.Variable(signal) # Variable in ordrer to modify at the next iteration
    nmem+=np.array(signal.shape).sum()
    hdata  = tf.constant(weights)
    nmem+=np.array(weights.shape).sum()

    #bidx = tf.constant(idx)
    #tf_bidx = tf.constant(idx//(netinfo["FSL"]['XIMAGE_SIZE']*netinfo["FSL"]['YIMAGE_SIZE']))

    # A MODIFIER HIDX
    #th,ph=hp.pix2ang(2048,realpix[hidx])
    #tmp_idx=hp.ang2pix(int(np.sqrt(netinfo["MAP"]['XIMAGE_SIZE']//12)),th,ph,nest=True)
    scale=(12*2048*2048)//netinfo["MAP"]['XIMAGE_SIZE']
    tmp_idx=realpix[hidx]//scale
    #del(th)
    #del(ph)
    tf_hidx = tf.constant(hidx.astype('int32'))
    nmem+=np.array(hidx.shape).sum()
    tf_idx_q = tf.constant((2*tmp_idx).astype('int32'))
    nmem+=np.array(tmp_idx.shape).sum()
    tf_idx_u = tf.constant((2*tmp_idx+1).astype('int32'))
    nmem+=np.array(tmp_idx.shape).sum()
    tf_weight_map = tf.constant(iiw_signal.astype("float32"))
    nmem+=np.array(iiw_signal.shape).sum()

    #modif neural network weight and biases
    for elt in netinfo:
      netinfo[elt]['dconv_weights'] = mynetwork[elt]['dconv_weights']
      netinfo[elt]['dconv_biases']  = mynetwork[elt]['dconv_biases']

    if mynetwork['MAP']['LEARNPARAM']==False and mynetwork['MAP']['LEARNFL']==False and mynetwork['MAP']['LEARNCONV']==False and mynetwork['MAP']['NFC']==False:
      themap = tf.Variable(np.zeros([netinfo["MAP"]['XIMAGE_SIZE']*2],dtype='float32'),name='MAP')
    else:
      themap = tf.reshape(model(params["MAP"],netinfo["MAP"],mynetwork["MAP"]["rank"]),[netinfo["MAP"]['XIMAGE_SIZE']*2])
      #themap = mynetwork['MAP']['RAP']*tf.reshape(model(params["MAP"],netinfo["MAP"],mynetwork["MAP"]["rank"]),[netinfo["MAP"]['XIMAGE_SIZE']*2])
    
    ampfsl= mynetwork['FSL']['RAP']*np.ones([netinfo["FSL"]['BATCHSZ'],netinfo["FSL"]['XIMAGE_SIZE']*netinfo["FSL"]['YIMAGE_SIZE']],dtype='float32')
    for i in range(netinfo["FSL"]['BATCHSZ']):
      ampfsl[i,:]*=coef[i]
  
    
    tf_ampfsl=tf.constant(ampfsl.reshape(netinfo["FSL"]['BATCHSZ']*netinfo["FSL"]['XIMAGE_SIZE']*netinfo["FSL"]['YIMAGE_SIZE']))
    nmem+=np.array(ampfsl.shape).sum()
    
    if mynetwork['FSL']['LEARNPARAM']==False and mynetwork['FSL']['LEARNFL']==False and mynetwork['FSL']['LEARNCONV']==False and mynetwork['FSL']['NFC']==False:
      thecorrfsl = tf.Variable(np.zeros([netinfo["FSL"]['BATCHSZ']*netinfo["FSL"]['XIMAGE_SIZE']*netinfo["FSL"]['YIMAGE_SIZE']],dtype='float32'),name='FSL')
      corrfsl =  tf_ampfsl*thecorrfsl
    else:
      corrfsl = tf_ampfsl*tf.reshape(model(params["FSL"],netinfo["FSL"],mynetwork["FSL"]['rank']),[netinfo["FSL"]['BATCHSZ']*netinfo["FSL"]['XIMAGE_SIZE']*netinfo["FSL"]['YIMAGE_SIZE']])
      
    toifsl = tf.gather(corrfsl,idx)
    
    imfsl = tf.math.unsorted_segment_sum(hdata*toifsl, tf_hidx, realpix.shape[0])
    deltag= (toifsl - tf.gather(tf_weight_map*imfsl,tf_hidx))
    
    ampmap = tf.constant([1.0]) 

    vsignal_q = ampmap*tf.gather(themap,2*tmp_idx) #tf_idx_q)
    vsignal_u = ampmap*tf.gather(themap,2*tmp_idx+1) #tf_idx_u)

    outmap = ampmap*themap
    
    out_signal = {}
    out_signal[0] =  toifsl + tf_OUT_CO*vsignal_q + tf_OUT_SI*vsignal_u
    out_signal[1] =  tf_OUT_CO*vsignal_q + tf_OUT_SI*vsignal_u
    out_signal[2] =  -tf_OUT_SI*vsignal_q + tf_OUT_CO*vsignal_u
    
    if mynetwork["MAP"]['rank']==0:
      sys.stderr.write('F MEM %.2fGB AVAIL %.2fGB\n'%(ps.virtual_memory().used/1024/1024/1024,ps.virtual_memory().available/1024/1024/1024))
        
    #pred = wu*(tf_CO*vsignal_q + tf_SI*vsignal_u)+wq*(-tf_SI*vsignal_q + tf_CO*vsignal_u)+deltag
    pred= vsignal_q*tf_CO + vsignal_u*tf_SI+deltag

    # *******************************************************************

    loss_type = mynetwork['CNN_LOSS']

    if loss_type=='RMSE':
      loss=tf.reduce_sum(hdata*tf.square(pred-data))  #PRED COMMON FACTOR!  
    elif loss_type=='EXP':
      loss = tf.reduce_sum(1-tf.exp(-(hdata*((pred-data)**2))))
    elif loss_type=='LOGCOSH':
      loss=tf.reduce_sum(tf.math.log(tf.math.cosh(hdata*(pred-data)))) #PRED COMMON FACTOR
    elif loss_type=='HUBER':
      loss = tf.losses.huber_loss(hdata*data, hdata*pred,delta=mynetwork['CNN_HUBER_LOSS_DELTA'])
    else:
      loss=tf.reduce_sum(hdata*tf.square(pred-data)) #PRED COMMON FACTOR!  
      if mynetwork['MAP']['rank'] == 0:
        print('WARNING: Unrecognized loss selection. Using default RMSE loss.')

    loss=loss #+1E6*loss3
    
    numbatch = tf.Variable(0, dtype=tf.float32)

    # Decay once per epoch, using an exponential schedule starting at 0.01.
    learning_rate = tf.train.exponential_decay(
        mynetwork['MAP']['LEARNING_RATE'],         # Base learning rate.
        numbatch,                           # Current index into the dataset.
        10,                                 # Decay step.
        mynetwork['MAP']['DECAY_RATE'],            # Decay rate.
        staircase=True)

    # Use simple momentum for the optimization.
    opti=tf.compat.v1.train.AdamOptimizer(learning_rate,0.9)
    optimizer = opti.minimize(loss,global_step=numbatch)
    #define igrad
    igrad = {}
    gradient={}
    assign = {}

    nnvar = 0
    
    if mynetwork['MAP']['LEARNPARAM']==False and mynetwork['MAP']['LEARNFL']==False and mynetwork['MAP']['LEARNCONV']==False and mynetwork['MAP']['NFC']==False:
      gradient[nnvar] = opti.compute_gradients(loss,var_list=[themap])[0]
      igrad[nnvar] = tf.placeholder(tf.float32,shape=(netinfo["MAP"]['XIMAGE_SIZE']*2))
      nnvar+=1
    if mynetwork['FSL']['LEARNPARAM']==False and mynetwork['FSL']['LEARNFL']==False and mynetwork['FSL']['LEARNCONV']==False and mynetwork['FSL']['NFC']==False:
      gradient[nnvar] = opti.compute_gradients(loss,var_list=[thecorrfsl])[0]
      igrad[nnvar] = tf.placeholder(tf.float32,shape=(netinfo["FSL"]['BATCHSZ']*netinfo["FSL"]['XIMAGE_SIZE']*netinfo["FSL"]['YIMAGE_SIZE']))
      nnvar+=1
          
    for elt in netinfo:
      if mynetwork[elt]['LEARNPARAM']==True:
        gradient[nnvar] = opti.compute_gradients(loss,var_list=[params[elt]])[0]
        igrad[nnvar] = tf.placeholder(tf.float32,shape=(mynetwork[elt]['BATCHSZ'],mynetwork[elt]['NTOTHIDDEN']))
        nnvar+=1

      if mynetwork[elt]['LEARNFL']==True:
        gradient[nnvar] = opti.compute_gradients(loss,var_list=[mynetwork[elt]['fc2_weights']])[0]
        gradient[nnvar+1] = opti.compute_gradients(loss,var_list=[mynetwork[elt]['fc2_biases']])[0]

        igrad[nnvar] = tf.placeholder(tf.float32,shape=(mynetwork[elt]['NTOTHIDDEN'],mynetwork[elt]['XIMAGE_SIZE'] // (mynetwork[elt]['XSCALE']**mynetwork[elt]['NDCONV'])*mynetwork[elt]['YIMAGE_SIZE']// (mynetwork[elt]['YSCALE']**(mynetwork[elt]['NDCONV']))*mynetwork[elt]['NUM_DCONV_CHAN'][0]))
        igrad[nnvar+1] = tf.placeholder(tf.float32,shape=(mynetwork[elt]['XIMAGE_SIZE'] // (mynetwork[elt]['XSCALE']**mynetwork[elt]['NDCONV'])*mynetwork[elt]['YIMAGE_SIZE']// (mynetwork[elt]['YSCALE']**(mynetwork[elt]['NDCONV']))*mynetwork[elt]['NUM_DCONV_CHAN'][0]))
        nnvar+=2

      if mynetwork[elt]['LEARNCONV'] == True:
        for i in range(netinfo[elt]['NDCONV']):
          igrad[nnvar]    = tf.compat.v1.placeholder(tf.float32,shape=(netinfo[elt]['XKERNELSZ'],netinfo[elt]['YKERNELSZ'],
                                                                      mynetwork[elt]['NUM_DCONV_CHAN'][i+1], mynetwork[elt]['NUM_DCONV_CHAN'][i]))
          igrad[nnvar+1]  = tf.compat.v1.placeholder(tf.float32,shape=(mynetwork[elt]['NUM_DCONV_CHAN'][i+1]))

          gradient[nnvar] = opti.compute_gradients(loss,var_list=[netinfo[elt]['dconv_weights'][i]])[0]
          gradient[nnvar+1] = opti.compute_gradients(loss,var_list=[netinfo[elt]['dconv_biases'][i]])[0]

          nnvar+=2

      if mynetwork[elt]['NFC']==True:
        igrad[nnvar] = tf.compat.v1.placeholder(tf.float32,shape=(mynetwork[elt]['BATCHSZ'],mynetwork[elt]['XIMAGE_SIZE'] // (mynetwork[elt]['XSCALE']**mynetwork[elt]['NDCONV'])*mynetwork[elt]['YIMAGE_SIZE']// (mynetwork[elt]['YSCALE']**(mynetwork[elt]['NDCONV']))*mynetwork[elt]['NUM_DCONV_CHAN'][0]))
        gradient[nnvar] = opti.compute_gradients(loss,var_list= [params[elt]])[0]
        nnvar+=1
  
    if mynetwork['MAP']['rank']==0:
      for i in range(nnvar):
        print(gradient[i][0],gradient[i][1])
    tgradient = [(igrad[i],gradient[i][1]) for i in range(nnvar)]
    apply_grad = opti.apply_gradients(tgradient,global_step=numbatch)

    if mynetwork['MAP']['rank']==0:
      sys.stderr.write('Initialized Network\n')

    #=========================================================================
    # TO BE CONFIGURED !!! TO LIMIT THE THREADS USAGE
    config = tf.ConfigProto()
    config.intra_op_parallelism_threads = 2
    config.inter_op_parallelism_threads = 2
    sess=tf.Session()

    tf.global_variables_initializer().run(session=sess)
    
    if mynetwork["MAP"]['rank']==0:
      sys.stderr.write('G MEM %.2fGB AVAIL %.2fGB\n'%(ps.virtual_memory().used/1024/1024/1024,ps.virtual_memory().available/1024/1024/1024))
        
    del(iiw_signal)
    del(tmp_idx)
    del(CO)
    del(SI)
    del(OUT_CO)
    del(OUT_SI)
    del(ampfsl)
    del(idx)
    
    if mynetwork["MAP"]['rank']==0:
      sys.stderr.write('H MEM %.2fGB AVAIL %.2fGB\n'%(ps.virtual_memory().used/1024/1024/1024,ps.virtual_memory().available/1024/1024/1024))
        
    if mynetwork["MAP"]['rank']==0:
      sys.stderr.write('Initialized Variables\n')

    myrun={}
    myrun['loss']=loss
    myrun['regul']=coefregul
    myrun['optimizer']=optimizer
    myrun['sess']=sess
    #myrun['logits']=themap
    myrun['learning_rate']=learning_rate
    myrun['params']=params
    #if mynetwork['DOFULLYCONNECTED']:
    #    myrun['flw']=netinfo['fc2_weights']
    #    myrun['flb']=netinfo['fc2_biases']
    #myrun['ISFL']=netinfo['DOFULLYCONNECTED']
    #myrun['dconvw']=dconv_weights
    #myrun['dconvb']=dconv_biases
    myrun['gradient']=gradient
    myrun['igrad']=igrad
    myrun['vsignal']=out_signal
    myrun['assign']=assign
    myrun['apply_grad']=apply_grad
    myrun['rank']=mynetwork['MAP']['rank']
    #myrun['RAPDATA']=mynetwork['MAP']['RAPDATA']
    myrun['FromFile']=mynetwork['MAP']['FromFile']
    myrun['nbolo']=netinfo['FSL']['nbolo']
    myrun['ngradient']=nnvar
    #myrun['SAVE_CNN_WEIGHTS']=mynetwork['MAP']['SAVE_CNN_WEIGHTS']
    #myrun['CNN_WEIGHTS']=mynetwork['MAP']['CNN_WEIGHTS']
    #myrun['DOFULLYCONNECTED']=mynetwork['DOFULLYCONNECTED']
    #myrun['NDCONV']=mynetwork['NDCONV']
    myrun['DOTRANSFERLEARNING']=mynetwork['MAP']['DOTRANSFERLEARNING']
    myrun['TMPNAME']=mynetwork['MAP']['TMPNAME']
    myrun["myMap"]=outmap
    myrun["corrfsl"]=corrfsl
    myrun["data"]= data  
    myrun["ampmap"]= tf.square(ampmap)
    #myrun["loss3"]= loss3
    #myrun["coef_fsl"]=coef[0:netinfo["FSL"]['BATCHSZ']]
    myrun['fsl_sz']=[netinfo["FSL"]['BATCHSZ'],netinfo["FSL"]['XIMAGE_SIZE'],netinfo["FSL"]['YIMAGE_SIZE']]
    print('MEMORY USED #SAMPLES',nmem)
    print('MEMORY USED GB',4*nmem/1024/1024/1024)
    return(myrun)

# ---------------------------------------------------------------------------------------- 
def init_network_data(myrun,signal,weights,hidx,realpix): 
    
    ii_signal  = np.bincount(hidx,weights=weights*signal,minlength=realpix.shape[0])
    iiw_signal = np.bincount(hidx,weights=weights,minlength=realpix.shape[0])
    
    lidx=np.where(iiw_signal>0)[0]

    ii_signal[lidx]/=iiw_signal[lidx]

    signal=signal-ii_signal[hidx]
    
    myrun["data"].assign(signal.astype('float32'))

    del(ii_signal)
    del(iiw_signal)
    del(lidx)

    return(myrun)

# ---------------------------------------------------------------------------------------- 
def calc_grad(myrun):
    sess=myrun['sess']

    #if myrun["rank"]==0:
    #  corrfsl = (sess.run([myrun["corrfsl"]])[0]).reshape(4,256,256)
    #  np.save('corrfsl.npy',corrfsl)
    #  exit(0)
    vvv=sess.run([myrun['gradient']])[0]
    
    nnvar= myrun['ngradient']
    #print(nnvar,vvv[0])
    vgrad={}
    for i in range(nnvar):
      if 'dense_shape' in dir(vvv[i][0]):
        vgrad["%03d"%(i)]=np.bincount(vvv[i][0][1],weights=vvv[i][0][0],minlength=vvv[i][0][2][0]).astype('float32')
      else:
        vgrad["%03d"%(i)] = np.array(vvv[i])[0]
    return(vgrad)
# ----------------------------------------------------------------------------------------
def calc_opti(myrun):
    sess=myrun['sess']
    sess.run([myrun['optimizer']])
    vvv=sess.run([myrun['params']])[0]

    nnvar=0
    vgrad={}

    vgrad["%03d"%(nnvar)] = np.array(vvv)
    nnvar+=1

    if myrun['DOFULLYCONNECTED']==True:
        vw,vb=sess.run([myrun['flw'],myrun['flb']])
        vgrad["%03d"%(nnvar)] = np.array(vw)
        vgrad["%03d"%(nnvar+1)] = np.array(vb)
        nnvar+=2

    ndconv=len(myrun['dconvw'])
    for i in range(ndconv):
        vw,vb=sess.run([myrun['dconvw'][i],myrun['dconvb'][i]])
        vgrad["%03d"%(nnvar)] = np.array(vw)
        vgrad["%03d"%(nnvar+1)] = np.array(vb)
        nnvar+=2
    return(vgrad)
# ----------------------------------------------------------------------------------------
def get_loss(myrun):
    sess=myrun['sess']
    #sess=tf.Session()
    l,lr=sess.run([myrun['loss'],myrun['learning_rate']])
    res=np.zeros([4],dtype='float32')
    res[0]=l
    res[1]=myrun['regul']
    res[2]=lr
    res[3]=1.0 #should the amplitude of the map
    return(res)
# ---------------------------------------------------------------------------------------- 
def apply_grad(myrun,vgrad):
    sess=myrun['sess']
    feed_dict={}
    for i in range(len(vgrad)):
      #if i==0:
      #  print(vgrad["%03d"%(i)][0])
      #print(i,np.array(vgrad["%03d"%(i)]).shape)
      feed_dict[myrun['igrad'][i]]=vgrad["%03d"%(i)]

    sess.run(myrun['apply_grad'], feed_dict=feed_dict)
    
# ----------------------------------------------------------------------------------------
def apply_opti(myrun,vgrad):
    sess=myrun['sess']
    feed_dict={}
    for i in range(len(vgrad)):
        feed_dict[myrun['igrad'][i]]=vgrad["%03d"%(i)]
    sess.run(myrun['assign'], feed_dict=feed_dict)
# ----------------------------------------------------------------------------------------
def get_prediction(myrun,num):
  sess=myrun['sess']
  predictions=((sess.run([myrun['logits']])[0])).flatten()
  if ((myrun['rank']==0) and (myrun['SAVE_CNN_WEIGHTS']) and (not myrun['DOTRANSFERLEARNING'])):

    if myrun['DOFULLYCONNECTED']:
      fc2weights_out,fc2biases_out=sess.run([myrun['flw'],myrun['flb']])        

      np.save('%s_fc2_weights.npy'%(myrun['CNN_WEIGHTS']),np.array(fc2weights_out))
      np.save('%s_fc2_biases.npy'%(myrun['CNN_WEIGHTS']),np.array(fc2biases_out))

    for i in range(myrun['NDCONV']):
      weights_out,biases_out=sess.run([myrun['dconvw'][i],myrun['dconvb'][i]])

      #print(weights_out.shape)
      #print(biases_out.shape)

      np.save('%s_dofc_%d_weights_%d.npy'%(myrun['CNN_WEIGHTS'],myrun['DOFULLYCONNECTED'],i),np.array(weights_out))
      np.save('%s_dofc_%d_biases_%d.npy'%(myrun['CNN_WEIGHTS'],myrun['DOFULLYCONNECTED'],i),np.array(biases_out))
  return(predictions.astype('float32'))

# ---------------------------------------------------------------------------------------- 
def get_correction(myrun,num):
  sess=myrun['sess']

  predictions=(sess.run([myrun['vsignal'][num[0]]])[0])

  if myrun["rank"]==0 and num[0]==0:
    map,corrfsl = sess.run([myrun["myMap"],myrun["corrfsl"]])
    corrfsl=corrfsl.reshape(myrun["fsl_sz"][0],myrun["fsl_sz"][1],myrun["fsl_sz"][2])
    #for i in range(len(myrun["coef_fsl"])):
    #  corrfsl[i,:,:]*=myrun["coef_fsl"][i]
    np.save(myrun["TMPNAME"]+"myMap_out.npy",map)
    np.save(myrun["TMPNAME"]+"corrfsl.npy",corrfsl)

    
  return(predictions.astype('float32'))

# ----------------------------------------------------------------------------------------
def close_session(myrun):
  sess=myrun['sess']
  
  sess.close()
  
  return(sess)
# ----------------------------------------------------------------------------------------




