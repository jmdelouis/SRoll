import numpy as np
import sys
#import matplotlib.pyplot as plt
import healpy as hp
import time

import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # Only errors displayed
import tensorflow.compat.v1 as tf

tf.logging.set_verbosity(tf.logging.ERROR)  # Double check: only errors displayed

tf.disable_v2_behavior()

"""
#plt.switch_backend('agg') # Backend for plots without displaying them (OCCIGEN)

print("Num GPUs Available: ", len(tf.config.experimental.list_physical_devices('GPU')))
tf.debugging.set_log_device_placement(True)

gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
  try:
    # Currently, memory growth needs to be the same across GPUs
    for gpu in gpus:
      tf.config.experimental.set_memory_growth(gpu, True)
    logical_gpus = tf.config.experimental.list_logical_devices('GPU')
    print(len(gpus), "Physical GPUs,", len(logical_gpus), "Logical GPUs")
  except RuntimeError as e:
    # Memory growth must be set before GPUs have been initialized
    print(e)
"""

def plt_vec(data):
  #plt.plot(data)
  #plt.show()
  del(data)
 
def plt_histo(data,hdata):
  idx=((data*256)/(2*np.pi)).astype('int')
  idx[idx==256]=255

  hh=np.bincount(idx,weights=hdata,minlength=256)

  if hh.std()>200:
    print('VAR %.4g %.4g'%(hh.std(),hh.mean()))
    #plt.plot(data[hdata>0])
    #plt.show()
    #plt.plot(hh)
    #plt.show()
  del(data)
  del(hdata)
 
def lrelu(x, alpha=0.3):
  #return tf.nn.relu(x)
  return tf.maximum(x, tf.multiply(x, alpha))
 
def wdecod(parameter,netinfo,rank): 
  # convert into image/bias
  #hidden = tf.matmul(parameter, fc1_weights) + fc1_biases
  #if DORELU==True:
  #    hidden = lrelu(hidden)
  if rank==0:
    sys.stderr.write('NPAR PER IMAGE %d\n'%(int((parameter.get_shape().as_list())[1])))
  if netinfo['DOFULLYCONNECTED']==False:
    hidden = parameter
  else:
    hidden = tf.matmul(parameter, netinfo['fc2_weights']) + netinfo['fc2_biases']
    if netinfo['DORELU']==True:
      hidden = lrelu(hidden)
  #hidden = lrelu(hidden)
  #old version ## NTIME=int((parameter.get_shape().as_list())[0])
  NTIME=int((hidden.get_shape().as_list())[0])

  if rank==0:
    sys.stderr.write('NPAR SIZE %d %d\n'%(int((parameter.get_shape().as_list())[0]),int((parameter.get_shape().as_list())[1])))

  # old version ## relu = tf.reshape(hidden,[NTIME, netinfo['XIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']) , 1 , netinfo['NUM_DCONV_CHAN'][0]])
  relu = tf.reshape(hidden,[NTIME, netinfo['XIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']), netinfo['YIMAGE_SIZE'] , netinfo['NUM_DCONV_CHAN'][0]])

  ## deconvolve : First layer
  # DeConvolution kernel
  for i in range(netinfo['NDCONV']):
    conv = tf.nn.conv2d_transpose(relu,netinfo['dconv_weights'][i],strides=[1, netinfo['SCALE'], 1, 1],padding='SAME',
                                  output_shape=[NTIME,netinfo['XIMAGE_SIZE'] // (netinfo['SCALE']**(netinfo['NDCONV']-1-i)),
                                  netinfo['YIMAGE_SIZE'],netinfo['NUM_DCONV_CHAN'][i+1]],name='dconv_%d'%(i))

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

  relu=relu/10

  return relu
 
def model(parameter,netinfo,rank):

  modpar=wdecod(parameter,netinfo,rank)
  return modpar
 
def alloc_table_float32(ndata):
  mydata = np.zeros([ndata],dtype='float32')
  return mydata
 
def alloc_table_int32(ndata):
  mydata = np.zeros([ndata],dtype='int32')
  return mydata
 
def free_table(data):
  del(data)
 
def init_shape(nbolo,xsize,ysize,rank,
               ncomp,ndconv,kernelsz,scale,npar_fc,
               loss_choice,flag_normalize_data_cnn_loss,huber_loss_delta,
               num_epochs,eval_freq,learning_rate,decay_rate,batchsz,
               flag_do_transfer_learning,flag_do_fully_connected,flag_do_relu,
               flag_test_convergence,flag_no_net,
               flag_save_cnn_weights, path_cnn_weights,
               seed):
  
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

  # CNN loss function parameters
  mynetwork['CNN_LOSS'] = loss_choice   
  mynetwork['CNN_NORMALIZE_DATA_IN_LOSS'] = flag_normalize_data_cnn_loss
  mynetwork['CNN_HUBER_LOSS_DELTA'] = huber_loss_delta,
               

  # Training parameters
  mynetwork['NUM_EPOCHS'] = num_epochs
  EVAL_FREQUENCY = eval_freq
  mynetwork['LEARNING_RATE']= learning_rate
  mynetwork['DECAY_RATE']=decay_rate
  mynetwork['BATCHSZ'] = batchsz
  
  # Debug options
  TESTCONVERGENCE = flag_test_convergence
  nonet = flag_no_net
  
  # Trasnfer learning
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
      NUM_DCONV_CHAN[i]=int(DEEPNESS)

    NUM_DCONV_CHAN[0] =NHIDDEN

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
 
def init_shape_from_file(TFLEARN,bolotab,rank,tmpid,tmpname):
  INFO=np.load('%s%s_NET_INFO.npy'%(TFLEARN,tmpname))

  nbolo=int(INFO[0])
  xsize=int(INFO[1])
  ysize=int(INFO[2])
  NDCONV=int(INFO[3])
  DOFULLYCONNECTED=int(INFO[4])==1
  mynetwork={}
  mynetwork['FromFile']=True

  if DOFULLYCONNECTED==True:
    for ib in range(nbolo):
      mynetwork['fc2w_%d'%(ib)] = np.load('%s%s_%s_NET_fc2w.npy'%(TFLEARN,tmpname,bolotab[ib]))
      mynetwork['fc2b_%d'%(ib)] = np.load('%s%s_%s_NET_fc2b.npy'%(TFLEARN,tmpname,bolotab[ib]))


  for ib in range(nbolo):
    for i in range(NDCONV):
      mynetwork['lw%d'%(i+ib*NDCONV)]=np.load('%s%s_%s_NET_w_%d.npy'%(TFLEARN,tmpname,bolotab[ib],i))
      mynetwork['lb%d'%(i+ib*NDCONV)]=np.load('%s%s_%s_NET_b_%d.npy'%(TFLEARN,tmpname,bolotab[ib],i))

  mynetwork['XIMAGE_SIZE']=xsize
  mynetwork['YIMAGE_SIZE']=ysize
  mynetwork['nbolo']=nbolo
  mynetwork['rank']=rank
  mynetwork['NCOMP']=int(INFO[5])
  mynetwork['SCALE']=int(INFO[7])
  mynetwork['NDCONV']=NDCONV
  XNHIDDEN=mynetwork['XIMAGE_SIZE'] // (mynetwork['SCALE']**mynetwork['NDCONV'])
  YNHIDDEN=mynetwork['YIMAGE_SIZE']
  NHIDDEN=mynetwork['NCOMP']
  DEEPNESS=int(INFO[6])
  EVAL_FREQUENCY = 100
  mynetwork['KERNELSZ'] = int(INFO[8])

  mynetwork['NUM_EPOCHS']= 2000
  SEED = 1234
  mynetwork['DORELU'] = True
  mynetwork['LEARNING_RATE']= 0.03
  mynetwork['DECAY_RATE']=0.995
  mynetwork['DOFULLYCONNECTED'] = True
  mynetwork['BATCHSZ'] = 1

  TESTCONVERGENCE = False
  
  NUM_DCONV_CHAN = {}

  if netinfo['DOFULLYCONNECTED']==False:
        NTOTHIDDEN=NHIDDEN*(XNHIDDEN*YNHIDDEN)
        for i in range(netinfo['NDCONV']):
            #NUM_DCONV_CHAN[i]=int(DEEPNESS*2**i)
            NUM_DCONV_CHAN[i]=int(DEEPNESS)
        NUM_DCONV_CHAN[0] =NHIDDEN
  else:
        NTOTHIDDEN=NHIDDEN
        for i in range(netinfo['NDCONV']):
            #NUM_DCONV_CHAN[i]=DEEPNESS*2**i
            NUM_DCONV_CHAN[i]=DEEPNESS
  NUM_DCONV_CHAN[netinfo['NDCONV']]=1 #DIFFERENT!

  netinfo['NUM_DCONV_CHAN']=NUM_DCONV_CHAN
  netinfo['NTOTHIDDEN']=NTOTHIDDEN

  if mynetwork['rank']==0:
    sys.stderr.write('init shape done\n')
  return(mynetwork)
 
def alloc_param(mynetwork):
  param=np.zeros([mynetwork['BATCHSZ']*mynetwork['NTOTHIDDEN']],dtype='float32')
  return(param)
 
def get_param(myrun):
  sess=myrun['sess']
  params=sess.run([myrun['params']])[0]
  return(params)
 
def alloc_convw(mynetwork):
  outobj={}
  for i in range(mynetwork['NDCONV']):
    outobj["%03d"%(i)] = np.zeros([mynetwork['KERNELSZ'] *mynetwork['NUM_DCONV_CHAN'][i+1]*mynetwork['NUM_DCONV_CHAN'][i]],dtype='float32')
  return(outobj)
 
def alloc_convb(mynetwork):
  outobj={}
  for i in range(mynetwork['NDCONV']):
    outobj["%03d"%(i)] = np.zeros([mynetwork['NUM_DCONV_CHAN'][i+1]],dtype='float32')
  return(outobj)
 
def get_convw(myrun):
  sess=myrun['sess']
  outobj={}
  for i in range(mynetwork['NDCONV']):
    outobj["%03d"%(i)] = sess.run([myrun['dconvw']])[0]
  return(outobj)
 
def get_convb(myrun):
  sess=myrun['sess']
  outobj={}
  for i in range(mynetwork['NDCONV']):
    outobj["%03d"%(i)] = sess.run([myrun['dconvb']])[0]
  return(outobj)
 
def alloc_flw(mynetwork):
  flw={}
  ib=0
  flw['%03d'%(ib)]=np.zeros([mynetwork['NTOTHIDDEN']* mynetwork['XIMAGE_SIZE'] // (mynetwork['SCALE']**mynetwork['NDCONV']) * mynetwork['YIMAGE_SIZE'] * mynetwork['NUM_DCONV_CHAN'][0]],dtype='float32')
  return(flw)
 
def alloc_flb(mynetwork):
  flb={}
  ib=0
  flb['%03d'%(ib)]=np.zeros([mynetwork['XIMAGE_SIZE'] // (mynetwork['SCALE']**mynetwork['NDCONV']) * mynetwork['YIMAGE_SIZE'] * mynetwork['NUM_DCONV_CHAN'][0]],dtype='float32')
  return(flb)
 
def get_flw(myrun):
  sess=myrun['sess']
  flw= sess.run([myrun['fc2_weights']])[0]
  return(flw)
 
def get_flb(myrun):
  sess=myrun['sess']
  flb= sess.run([myrun['fc2_biases']])[0]
  return(flb)
 
def init_network(mynetwork,signal,weights,TCO1,TSI1,MAT0,MAT1,MAT2,hidx,idx,realpix,in_param,in_convw,in_convb,in_flw,in_flb):
    
    mapq=np.load('/home1/datawork/tfoulqui/data/MAP_JMD_2048/DUSTMODEL_857_Q.npy').astype('float32')
    mapu=np.load('/home1/datawork/tfoulqui/data/MAP_JMD_2048/DUSTMODEL_857_U.npy').astype('float32')

    nonet = False
    RMAX = 1000000
    baserandom = np.random.rand(RMAX).astype('float32')-0.5
    MAXRAND = baserandom.shape[0]
    netinfo={}

    for i in mynetwork:
      netinfo[i]=mynetwork[i]
    NUM_DCONV_CHAN=mynetwork['NUM_DCONV_CHAN']

    if mynetwork['rank']==0:
      sys.stderr.write('VARIANCE CNN INPUT %f %f\n'%((signal*mynetwork['RAPDATA']).std(),(weights/mynetwork['RAPDATA']).std()))

    NPIX=realpix.shape[0]

    coefregul=weights.sum()

    #params = tf.Variable(in_param.reshape(mynetwork['BATCHSZ'],mynetwork['NTOTHIDDEN'],1))
    if nonet:
        params = tf.Variable(0.1*baserandom[np.arange(netinfo['BATCHSZ']*netinfo['XIMAGE_SIZE']*netinfo['YIMAGE_SIZE'],dtype='int')%MAXRAND].reshape([netinfo['BATCHSZ']*netinfo['XIMAGE_SIZE']*netinfo['YIMAGE_SIZE']]))

    else:
        params = tf.Variable(0.1*baserandom[np.arange(netinfo['BATCHSZ']*netinfo['NTOTHIDDEN'],dtype='int')%MAXRAND].reshape(netinfo['BATCHSZ'],netinfo['NTOTHIDDEN']))


    ii_signal=np.bincount(hidx,weights=weights*signal,minlength=realpix.shape[0])
    iiw_signal=np.bincount(hidx,weights=weights,minlength=realpix.shape[0])
    iico1 = np.bincount(hidx,weights=weights*TCO1,minlength=realpix.shape[0])
    iisi1 = np.bincount(hidx,weights=weights*TSI1,minlength=realpix.shape[0])
    lidx=np.where(iiw_signal>0)[0]

    ii_signal[lidx]/=iiw_signal[lidx]
    iico1[lidx]/=iiw_signal[lidx]
    iisi1[lidx]/=iiw_signal[lidx]


    signal=signal-ii_signal[hidx]


    tf_OUT_CO = tf.constant(TCO1.astype('float32'))
    tf_OUT_SI = tf.constant(TSI1.astype('float32'))

    TCO1 = (TCO1 - iico1[hidx]).astype('float32')
    TSI1 = (TSI1 - iisi1[hidx]).astype('float32')

    data   = tf.constant(signal.astype('float32'))

    hdata  = tf.constant(weights)

    tf_CO  = tf.constant(TCO1)

    tf_SI  = tf.constant(TSI1)

    # A MODIFIER HIDX
    th,ph=hp.pix2ang(2048,realpix[hidx])
    idx=hp.ang2pix(int(np.sqrt(netinfo['XIMAGE_SIZE']//12)),th,ph,nest=True)
    tf_hidx = tf.constant(hidx.astype('int32'))

    tf_idx_q = tf.constant((2*idx).astype('int32'))
    tf_idx_u = tf.constant((2*idx+1).astype('int32'))

    if mynetwork['DOFULLYCONNECTED']:
      if mynetwork['DOTRANSFERLEARNING']:
        netinfo['fc2_weights'] = tf.Variable(np.load('%s_fc2_weights.npy'%(mynetwork['CNN_WEIGHTS'])))
        netinfo['fc2_biases'] = tf.Variable(np.load('%s_fc2_biases.npy'%(mynetwork['CNN_WEIGHTS'])))
      else:
        #netinfo['fc2_weights'] = tf.Variable(in_flw[0].reshape(NTOTHIDDEN, netinfo['XIMAGE_SIZE']// (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['YIMAGE_SIZE'] * netinfo['NUM_DCONV_CHAN'][0]))
        #netinfo['fc2_biases'] = tf.Variable(in_flb[0].reshape(netinfo['XIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['YIMAGE_SIZE'] * netinfo['NUM_DCONV_CHAN'][0]))
        netinfo['fc2_weights'] = tf.Variable((np.random.rand(netinfo['NTOTHIDDEN'], netinfo['XIMAGE_SIZE']// (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['YIMAGE_SIZE'] * netinfo['NUM_DCONV_CHAN'][0])-0.5).astype('float32'))
        netinfo['fc2_biases'] = tf.Variable((np.random.rand(netinfo['XIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['YIMAGE_SIZE'] * netinfo['NUM_DCONV_CHAN'][0])-0.5).astype('float32'))

    '''
    # BASE FOR NEXT FSL CNN IMPLEMENTATION
    if netinfo['DOFULLYCONNECTED']:
      if mynetwork['FromFile']==True:
        netinfo['fc2_weights'] = {}
        netinfo['fc2_biases'] = {}
        netinfo['fc2_weights'] = tf.constant(mynetwork['fc2w_%d'%(ib)])
        netinfo['fc2_biases'] = tf.constant(mynetwork['fc2b_%d'%(ib)])
      else:
        ib=0
        netinfo['fc2_weights_%d'%(ib)] = tf.Variable(in_flw['%03d'%(ib)].reshape(NTOTHIDDEN, netinfo['XIMAGE_SIZE']// (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['YIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['NUM_DCONV_CHAN'][0]))
        netinfo['fc2_biases_%d'%(ib)] = tf.Variable(in_flb['%03d'%(ib)].reshape(netinfo['XIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['YIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['NUM_DCONV_CHAN'][0]))
    '''

    #modif neural network weight and biases
    #Init
    dconv_weights={}
    dconv_biases={}

    for i in range(netinfo['NDCONV']):
        
        if mynetwork['DOTRANSFERLEARNING']:
          dconv_weights[i] = tf.Variable(np.load('%s_dofc_%d_weights_%d.npy'%(mynetwork['CNN_WEIGHTS'],mynetwork['DOFULLYCONNECTED'],i)))
          dconv_biases[i] = tf.Variable(np.load('%s_dofc_%d_biases_%d.npy'%(mynetwork['CNN_WEIGHTS'],mynetwork['DOFULLYCONNECTED'],i)))
        else:        
          #print('INIT WEIGHTS ')
          tmp=np.random.randn(mynetwork['KERNELSZ'], 1,NUM_DCONV_CHAN[i+1], NUM_DCONV_CHAN[i]).astype('float32')
          #Calcul weight et biase pour 1ere et derniere couche
          if i==0 or i== netinfo['NDCONV']-1:
              dconv_weights[i] = tf.Variable(tmp)
              tmp=np.zeros(NUM_DCONV_CHAN[i+1]).astype('float32')
              dconv_biases[i] = tf.Variable(tmp)
          #Sinon weight = weight(couche-1) and biase= biase(couche-1)
          else:
              dconv_weights[i] = dconv_weights[i-1]
              dconv_biases[i] = dconv_biases[i-1]


    netinfo['dconv_weights']=dconv_weights
    netinfo['dconv_biases']=dconv_biases

    if nonet:
        logits = params
    else:
        lpar = params
        logits = tf.reshape(model(lpar,netinfo,mynetwork['rank']),[netinfo['BATCHSZ']*netinfo['XIMAGE_SIZE']*netinfo['YIMAGE_SIZE']])


    #logits = tf.reshape(model(lpar,netinfo,mynetwork['rank']),[netinfo['BATCHSZ']*netinfo['XIMAGE_SIZE']*netinfo['YIMAGE_SIZE']])

    vsignal_q = tf.gather(logits,tf_idx_q)
    vsignal_u = tf.gather(logits,tf_idx_u)

    vmap_q = tf.gather(mapq,tf_idx_q)
    vmap_u = tf.gather(mapu,tf_idx_q)

    pred=tf_CO*vsignal_q + tf_SI*vsignal_u

    out_signal = tf_OUT_CO*vsignal_q + tf_OUT_SI*vsignal_u
    out_sigmap = tf_OUT_CO*vmap_q + tf_OUT_SI*vmap_u

    if mynetwork['CNN_NORMALIZE_DATA_IN_LOSS']:
      data_std = tf.sqrt(tf.reduce_sum(hdata*((data)**2)))
      data = data/data_std
      pred = pred/data_std
 
    loss_type = mynetwork['CNN_LOSS']

    if loss_type=='RMSE':
      loss=tf.reduce_sum(hdata*pred*(pred-2*data)) #PRED COMMON FACTOR!  
    elif loss_type=='EXP':
      loss = tf.reduce_sum(1-tf.exp(-(hdata*((pred-data)**2))))
    elif loss_type=='LOGCOSH':
      loss=tf.reduce_sum(tf.math.log(tf.math.cosh(hdata*(pred-data)))) #PRED COMMON FACTOR
    elif loss_type=='HUBER':
      loss = tf.losses.huber_loss(hdata*data, hdata*pred,delta=mynetwork['CNN_HUBER_LOSS_DELTA'])
    else:
      loss=tf.reduce_sum(hdata*pred*(pred-2*data)) #PRED COMMON FACTOR!  
      print('WARNING: Unrecognized loss selection. Using default RMSE loss.')

    numbatch = tf.Variable(0, dtype=tf.float32)

    # Decay once per epoch, using an exponential schedule starting at 0.01.
    learning_rate = tf.train.exponential_decay(
        mynetwork['LEARNING_RATE'],         # Base learning rate.
        numbatch,                           # Current index into the dataset.
        10,                                 # Decay step.
        mynetwork['DECAY_RATE'],            # Decay rate.
        staircase=True)

    # Use simple momentum for the optimization.
    opti=tf.compat.v1.train.AdamOptimizer(learning_rate,0.9)
    optimizer = opti.minimize(loss,global_step=numbatch)
    #define igrad
    igrad = {}
    gradient={}
    assign = {}

    nnvar = 0

     #igrad[0]    = tf.compat.v1.placeholder(tf.float32,shape=(netinfo['BATCHSZ'],netinfo['NTOTHIDDEN'],1))
    if nonet:
        igrad[nnvar]    = tf.placeholder(tf.float32,shape=(netinfo['BATCHSZ']*netinfo['XIMAGE_SIZE']*netinfo['YIMAGE_SIZE']))
    else:
        igrad[nnvar]    = tf.placeholder(tf.float32,shape=(netinfo['BATCHSZ'],netinfo['NTOTHIDDEN']))
  

    gradient[nnvar] = opti.compute_gradients(loss,var_list=[params])[0]
    assign[nnvar]   = params.assign(igrad[nnvar])

    nnvar+=1

    if mynetwork['DOTRANSFERLEARNING']==False:
  
      if mynetwork['DOFULLYCONNECTED']:

        igrad[nnvar]    = tf.compat.v1.placeholder(tf.float32,shape=(netinfo['NTOTHIDDEN'], netinfo['XIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV'])* netinfo['YIMAGE_SIZE'] * netinfo['NUM_DCONV_CHAN'][0]))
        igrad[nnvar+1]    = tf.compat.v1.placeholder(tf.float32,shape=(netinfo['XIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['YIMAGE_SIZE'] * netinfo['NUM_DCONV_CHAN'][0]))
        gradient[nnvar] = opti.compute_gradients(loss,var_list=[netinfo['fc2_weights']])[0]
        gradient[nnvar+1] = opti.compute_gradients(loss,var_list=[netinfo['fc2_biases']])[0]

        assign[nnvar] = netinfo['fc2_weights'].assign(igrad[nnvar])
        assign[nnvar+1] = netinfo['fc2_biases'].assign(igrad[nnvar+1])
        nnvar+=2

      for i in range(netinfo['NDCONV']):
        igrad[nnvar]    = tf.compat.v1.placeholder(tf.float32,shape=(netinfo['KERNELSZ'], 1, NUM_DCONV_CHAN[i+1], NUM_DCONV_CHAN[i]))
        igrad[nnvar+1]  = tf.compat.v1.placeholder(tf.float32,shape=(NUM_DCONV_CHAN[i+1]))

        gradient[nnvar] = opti.compute_gradients(loss,var_list=[dconv_weights[i]])[0]
        gradient[nnvar+1] = opti.compute_gradients(loss,var_list=[dconv_biases[i]])[0]

        assign[nnvar]   = dconv_weights[i].assign(igrad[nnvar])
        assign[nnvar+1] = dconv_biases[i].assign(igrad[nnvar+1])
        nnvar+=2

    tgradient = [(igrad[i],gradient[i][1]) for i in range(nnvar)]
    apply_grad = opti.apply_gradients(tgradient,global_step=numbatch)

    if mynetwork['rank']==0:
      sys.stderr.write('Initialized Network\n')

    sess=tf.Session()

    tf.global_variables_initializer().run(session=sess)
 
 
     
    #del
    del(signal)
    del(weights)
    del(TCO1)
    del(TSI1)
    del(MAT1)
    del(MAT2)
    del(MAT0)
    del(hidx)
    del(idx)
    #del(nhidx)
    del(realpix)
    #del(pidx)

 
 
 
    if mynetwork['rank']==0:
      sys.stderr.write('Initialized Variables\n')

    myrun={}
    myrun['loss']=loss
    myrun['regul']=coefregul
    myrun['optimizer']=optimizer
    myrun['sess']=sess
    myrun['logits']=logits
    myrun['learning_rate']=learning_rate
    myrun['params']=params
    if mynetwork['DOFULLYCONNECTED']:
        myrun['flw']=netinfo['fc2_weights']
        myrun['flb']=netinfo['fc2_biases']
    myrun['ISFL']=netinfo['DOFULLYCONNECTED']
    myrun['dconvw']=dconv_weights
    myrun['dconvb']=dconv_biases
    myrun['gradient']=gradient
    myrun['igrad']=igrad
    myrun['vsignal']=out_signal
    myrun['assign']=assign
    myrun['apply_grad']=apply_grad
    myrun['rank']=mynetwork['rank']
    myrun['RAPDATA']=mynetwork['RAPDATA']
    myrun['FromFile']=mynetwork['FromFile']
    myrun['nbolo']=netinfo['nbolo']
    myrun['ngradient']=nnvar
    myrun['SAVE_CNN_WEIGHTS']=mynetwork['SAVE_CNN_WEIGHTS']
    myrun['CNN_WEIGHTS']=mynetwork['CNN_WEIGHTS']
    myrun['DOFULLYCONNECTED']=mynetwork['DOFULLYCONNECTED']
    myrun['NDCONV']=mynetwork['NDCONV']
    myrun['DOTRANSFERLEARNING']=mynetwork['DOTRANSFERLEARNING']
    
    

    return(myrun)
 
def calc_grad(myrun):
    sess=myrun['sess']
    sess.run([myrun['optimizer']])
    vvv=sess.run([myrun['params']])[0]

    nnvar=0
    vgrad={}

    vgrad["%03d"%(nnvar)] = np.array(vvv)
    nnvar+=1

    if myrun['DOTRANSFERLEARNING']==False:

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
 
def get_loss(myrun):
    sess=myrun['sess']
    #sess=tf.Session()
    l,lr=sess.run([myrun['loss'],myrun['learning_rate']])
    res=np.zeros([2],dtype='float32')
    res[0]=l/myrun['regul']
    res[1]=lr
    return(res)
 
def apply_grad(myrun,vgrad):
    sess=myrun['sess']
    feed_dict={}
    for i in range(len(vgrad)):
        feed_dict[myrun['igrad'][i]]=vgrad["%03d"%(i)]
    sess.run(myrun['apply_grad'], feed_dict=feed_dict)
 
def apply_opti(myrun,vgrad):
    sess=myrun['sess']
    feed_dict={}
    for i in range(len(vgrad)):
        feed_dict[myrun['igrad'][i]]=vgrad["%03d"%(i)]
    sess.run(myrun['assign'], feed_dict=feed_dict)
'''

def get_prediction(myrun,output_name):
    sess=myrun['sess']

    #Create output name for POLARMAP file
    output_name = output_name.split("/")
    id_name = len(output_name)-1
    output_name= output_name[id_name]

    polmaps=sess.run([myrun['logits']])[0]
    predictions=sess.run([myrun['out_signal']])[0]
    predmap=sess.run([myrun['out_sigmap']])[0]


    if myrun['rank']==0:
        np.save(str(output_name)+'_POLARMAP.npy',polmaps)
        np.save(str(output_name)+'_OUT.npy',predictions)
        np.save(str(output_name)+'_OUTMAP.npy',predmap)

    if TESTCONVERGENCE == True:
        sys.stderr.write('Save the correction')
        np.save('predictions_theo_final_ndconv_7_fc.npy',predictions)


    sess.close()

    return(polmaps,predictions,predmap)
'''

def get_prediction(myrun):
  sess=myrun['sess']
  predictions=((sess.run([myrun['logits']])[0])/myrun['RAPDATA']).flatten()

  if ((myrun['rank']==0) and (myrun['SAVE_CNN_WEIGHTS']) and (not myrun['DOTRANSFERLEARNING'])):

    if myrun['DOFULLYCONNECTED']:
      fc2weights_out,fc2biases_out=sess.run([myrun['flw'],myrun['flb']])        

      np.save('%s_fc2_weights.npy'%(myrun['CNN_WEIGHTS']),np.array(fc2weights_out))
      np.save('%s_fc2_biases.npy'%(myrun['CNN_WEIGHTS']),np.array(fc2biases_out))

    for i in range(myrun['NDCONV']):
      weights_out,biases_out=sess.run([myrun['dconvw'][i],myrun['dconvb'][i]])

      print(weights_out.shape)
      print(biases_out.shape)

      np.save('%s_dofc_%d_weights_%d.npy'%(myrun['CNN_WEIGHTS'],myrun['DOFULLYCONNECTED'],i),np.array(weights_out))
      np.save('%s_dofc_%d_biases_%d.npy'%(myrun['CNN_WEIGHTS'],myrun['DOFULLYCONNECTED'],i),np.array(biases_out))

  return(predictions.astype('float32'))

 
def get_correction(myrun):
  sess=myrun['sess']
  predictions=(sess.run([myrun['vsignal']])[0])/myrun['RAPDATA']
  return(predictions.astype('float32'))
 
def close_session(myrun):
  sess=myrun['sess']
  sess.close()
  return(sess)
