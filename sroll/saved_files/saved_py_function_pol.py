import numpy as np
import sys
#import matplotlib.pyplot as plt
import healpy as hp
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()
import time



def plt_vec(data):
  plt.plot(data)
  plt.show()
  del(data)

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
  NTIME=int((parameter.get_shape().as_list())[0])
  if rank==0:
    sys.stderr.write('NPAR SIZE %d %d\n'%(int((parameter.get_shape().as_list())[0]),int((parameter.get_shape().as_list())[1])))

  relu = tf.reshape(hidden,[NTIME, netinfo['XIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']) , 1 , netinfo['NUM_DCONV_CHAN'][0]])
  #
  ## deconvolve : First layer
  # DeConvolution kernel
  for i in range(netinfo['NDCONV']):
    conv = tf.nn.conv2d_transpose(relu,netinfo['dconv_weights'][i],strides=[1, netinfo['SCALE'], 1, 1],padding='SAME',
                                  output_shape=[NTIME,netinfo['XIMAGE_SIZE'] // (netinfo['SCALE']**(netinfo['NDCONV']-1-i)),1,
                                                netinfo['NUM_DCONV_CHAN'][i+1]],name='dconv_%d'%(i))
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
    
    
def init_shape_from_file(TFLEARN,bolotab,rank,tmpid,tmpname):


    INFO=np.load('%s%s_NET_INFO.npy'%(TFLEARN,tmpname))
    
    nbolo=int(INFO[0]) 
    xsize=int(INFO[1])
    ysize=int(INFO[2])
    NDCONV=int(INFO[3])
    DOFULLYCONNECTED=int(INFO[4])==1
    
    TESTCONVERGENCE = False

    mynetwork={}

    mynetwork['FromFile']=True

    if DOFULLYCONNECTED==True:
      for ib in range(nbolo):
        mynetwork['fc2w_%d'%(ib)] = np.load('%s%s_%s_NET_fc2w.npy'%(TFLEARN,tmpname,bolotab[ib]))
        mynetwork['fc2b_%d'%(ib)] = np.load('%s%s_%s_NET_fc2b.npy'%(TFLEARN,tmpname,bolotab[ib]))
        
        
    for ib in range(nbolo):
      for i in range(NDCONV):
        mynetwork['lw%d_%d'%(i,ib)]=np.load('%s%s_%s_NET_w_%d.npy'%(TFLEARN,tmpname,bolotab[ib],i))
        mynetwork['lb%d_%d'%(i,ib)]=np.load('%s%s_%s_NET_b_%d.npy'%(TFLEARN,tmpname,bolotab[ib],i))
    
    mynetwork['XIMAGE_SIZE']=xsize
    mynetwork['YIMAGE_SIZE']=ysize
    mynetwork['nbolo']=nbolo
    mynetwork['rank']=rank
    mynetwork['NCOMP']=int(INFO[5])
    mynetwork['SCALE']=int(INFO[7])
    mynetwork['NDCONV']=NDCONV
    XNHIDDEN=mynetwork['XIMAGE_SIZE'] // (mynetwork['SCALE']**mynetwork['NDCONV'])
    YNHIDDEN=mynetwork['YIMAGE_SIZE'] // (mynetwork['SCALE']**mynetwork['NDCONV'])
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
    
    NUM_DCONV_CHAN = {}
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
    mynetwork['RAPDATA']=float(INFO[9])


    if mynetwork['rank']==0:
      sys.stderr.write('init shape done\n')
    return(mynetwork)

def init_shape(nbolo,xsize,ysize,rank,test_param):
    
    mynetwork={}
    mynetwork['xsize']=xsize
    mynetwork['ysize']=ysize
    mynetwork['nbolo']=nbolo
    mynetwork['rank']=rank
    
  
    if mynetwork['rank']==0:
      sys.stderr.write('init shape done\n')
    return(mynetwork)

def init_network(mynetwork,signal,weights,TCO1,TSI1,MAT0,MAT1,MAT2,hidx,idx,realpix,baserandom):
    netinfo={}
    netinfo['XIMAGE_SIZE']=mynetwork['xsize']
    netinfo['YIMAGE_SIZE']=mynetwork['ysize']
    netinfo['BATCHSZ']=1
    
    MAXRAND=baserandom.shape[0]
    
    if mynetwork['rank']==0:
      sys.stderr.write('VARIANCE CNN INPUT %f %f\n'%(signal.std(),weights.std()))

    NPIX=realpix.shape[0]
    netinfo['NCOMP']=8
    netinfo['SCALE']=4
    netinfo['NDCONV']=4
    XNHIDDEN=netinfo['XIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV'])
    YNHIDDEN=1
    NHIDDEN=netinfo['NCOMP']
    DEEPNESS=netinfo['NCOMP']
    EVAL_FREQUENCY = 100
    KERNELSZ = 4

    NUM_EPOCHS1= 0
    NUM_EPOCHS= 200
    SEED = 1234
    netinfo['DORELU'] = True
    LEARNING_RATE= 0.01
    DECAY_RATE=0.995
    netinfo['DOFULLYCONNECTED'] = False
    
    NUM_DCONV_CHAN = {}

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


    if netinfo['DOFULLYCONNECTED']==False:
        NTOTHIDDEN=NHIDDEN*(XNHIDDEN*YNHIDDEN)
        for i in range(netinfo['NDCONV']):
            NUM_DCONV_CHAN[i]=int(DEEPNESS*2**i)
        NUM_DCONV_CHAN[0] =NHIDDEN

    else:
        NTOTHIDDEN=1
        for i in range(netinfo['NDCONV']):
            NUM_DCONV_CHAN[i]=DEEPNESS*2**i
    NUM_DCONV_CHAN[netinfo['NDCONV']]=2

    netinfo['NUM_DCONV_CHAN']=NUM_DCONV_CHAN
    netinfo['NTOTHIDDEN']=NTOTHIDDEN

    coefregul=weights.sum()

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
    del(signal)
    hdata  = tf.constant(weights)
    del(weights)
    tf_CO  = tf.constant(TCO1)
    del(TCO1)
    tf_SI  = tf.constant(TSI1)
    del(TSI1)
    del(MAT1)
    del(MAT2)
    del(MAT0)
    # A MODIFIER HIDX
    th,ph=hp.pix2ang(2048,realpix[hidx])
    idx=hp.ang2pix(int(np.sqrt(netinfo['XIMAGE_SIZE']//12)),th,ph,nest=True)
    tf_hidx = tf.constant(hidx.astype('int32'))
    del(hidx)
    tf_idx_q = tf.constant((2*idx).astype('int32'))
    tf_idx_u = tf.constant((2*idx+1).astype('int32'))
    del(idx)

    if netinfo['DOFULLYCONNECTED']:
        netinfo['fc2_weights'] = tf.Variable(0.1*baserandom[np.arange(NTOTHIDDEN* netinfo['XIMAGE_SIZE']* (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['YIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['NUM_DCONV_CHAN'][0],dtype='int')%MAXRAND].reshape(NTOTHIDDEN, netinfo['XIMAGE_SIZE']// (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['YIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['NUM_DCONV_CHAN'][0]))
        netinfo['fc2_biases'] = tf.Variable(0.1*baserandom[np.arange(netinfo['XIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['YIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['NUM_DCONV_CHAN'][0],dtype='int')%MAXRAND].reshape(netinfo['XIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['YIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['NUM_DCONV_CHAN'][0]))

    dconv_weights={}
    dconv_biases={}
    nw=0
    for i in range(netinfo['NDCONV']):
        dconv_weights[i] = tf.Variable(0.1*baserandom[np.arange(KERNELSZ*NUM_DCONV_CHAN[i+1]*NUM_DCONV_CHAN[i],dtype='int')%MAXRAND].reshape(KERNELSZ,1,NUM_DCONV_CHAN[i+1],NUM_DCONV_CHAN[i]))
        dconv_biases[i] = tf.Variable(0.1*baserandom[np.arange(NUM_DCONV_CHAN[i+1],dtype='int')%MAXRAND])

    del(baserandom)
    del(realpix)

    netinfo['dconv_weights']=dconv_weights
    netinfo['dconv_biases']=dconv_biases

    lpar = params
    logits = tf.reshape(model(lpar,netinfo,mynetwork['rank']),[netinfo['BATCHSZ']*netinfo['XIMAGE_SIZE']*netinfo['YIMAGE_SIZE']])

    vsignal_q = tf.gather(logits,tf_idx_q)
    vsignal_u = tf.gather(logits,tf_idx_u)
   
    residu = data - (tf_CO*vsignal_q + tf_SI*vsignal_u)

    out_signal = tf_OUT_CO*vsignal_q + tf_OUT_SI*vsignal_u

    loss=tf.reduce_sum(hdata*tf.square(residu))

    numbatch = tf.Variable(0, dtype=tf.float32)

    # Decay once per epoch, using an exponential schedule starting at 0.01.
    learning_rate = tf.train.exponential_decay(
        LEARNING_RATE,       # Base learning rate.
        numbatch,            # Current index into the dataset.
        10,                  # Decay step.
        DECAY_RATE,          # Decay rate.
        staircase=True)

    # Use simple momentum for the optimization.
    opti=tf.compat.v1.train.AdamOptimizer(learning_rate,0.9)

    optimizer = opti.minimize(loss,global_step=numbatch)

    igrad={}   
    igrad[0]    = tf.compat.v1.placeholder(tf.float32,shape=(netinfo['BATCHSZ'],netinfo['NTOTHIDDEN']))
    assign = {}
    assign[0]   = params.assign(igrad[0])

    nnvar=1
    for i in range(netinfo['NDCONV']):
        igrad[nnvar+i*2]    = tf.compat.v1.placeholder(tf.float32,shape=(KERNELSZ, 1, NUM_DCONV_CHAN[i+1], NUM_DCONV_CHAN[i]))
        igrad[nnvar+i*2+1]  = tf.compat.v1.placeholder(tf.float32,shape=(NUM_DCONV_CHAN[i+1]))
        assign[nnvar+i*2]   = dconv_weights[i].assign(igrad[nnvar+i*2])
        assign[nnvar+i*2+1] = dconv_biases[i].assign(igrad[nnvar+i*2+1])
    
    if mynetwork['rank']==0:
      sys.stderr.write('Initialized Network\n')

    sess=tf.Session()

    tf.global_variables_initializer().run(session=sess)

    if mynetwork['rank']==0:
      sys.stderr.write('Initialized Variables\n')

    myrun={}
    myrun['loss']=loss
    myrun['regul']=coefregul
    myrun['optimizer']=optimizer
    myrun['sess']=sess
    myrun['logits']=logits
    myrun['out_signal']=out_signal
    myrun['learning_rate']=learning_rate
    myrun['params']=params
    if netinfo['DOFULLYCONNECTED']:
        myrun['flw']=netinfo['fc2_weights']
        myrun['flb']=netinfo['fc2_biases']
    myrun['ISFL']=netinfo['DOFULLYCONNECTED']
    myrun['dconvw']=dconv_weights
    myrun['dconvb']=dconv_biases
    myrun['igrad']=igrad
    myrun['assign']=assign
    myrun['rank']=mynetwork['rank']
    return(myrun)

def calc_grad(myrun):
    sess=myrun['sess']
    sess.run([myrun['optimizer']])
    vvv=sess.run([myrun['params']])[0]
    vgrad={}
    vgrad["%03d"%(0)] = np.array(vvv)
    if myrun['ISFL']==True:
        vw,vb=sess.run([myrun['flw'],myrun['flb']])
        vgrad["%03d"%(1)] = np.array(vw)
        vgrad["%03d"%(2)] = np.array(vb)
        nnvar=3
    else:
        nnvar=1
    ndconv=len(myrun['dconvw'])
    for i in range(ndconv):
        vw,vb=sess.run([myrun['dconvw'][i],myrun['dconvb'][i]])
        vgrad["%03d"%(nnvar+i*2)] = np.array(vw)
        vgrad["%03d"%(nnvar+i*2+1)] = np.array(vb)
    return(vgrad)

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

def alloc_flw(netinfo):
  flw={}
  ib=0
  flw['%03d'%(ib)]=np.zeros([NTOTHIDDEN* netinfo['XIMAGE_SIZE']* (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['YIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['NUM_DCONV_CHAN'][0]],dtype='float32')
  return(flw)

def alloc_flb(mynetwork):
  flb={}
  ib=0
  flb['%03d'%(ib)]=np.zeros([netinfo['XIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['YIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['NUM_DCONV_CHAN'][0]],dtype='float32')
  return(flb)

def get_flw(myrun):
  sess=myrun['sess']
  flw= sess.run([myrun['fc2_weights']])[0]
  return(flw)

def get_flb(myrun):
  sess=myrun['sess']
  flb= sess.run([myrun['fc2_biases']])[0]
  return(flb)
 
 
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

    sess.run(myrun['assign'], feed_dict=feed_dict)

def get_prediction(myrun):
    sess=myrun['sess']
    predictions=sess.run([myrun['logits']])[0]
    if myrun['rank']==0:
        np.save('POLARMAP.npy',predictions)
    predictions=sess.run([myrun['out_signal']])[0]
    if TESTCONVERGENCE == True:
      sys.stderr.write('Save the correction')
      np.save('predictions_new.npy',predictions)
    sess.close()
    return(predictions)   

def get_correction(myrun):
  sess=myrun['sess']
  predictions=(sess.run([myrun['logits']])[0])/myrun['RAPDATA']
  return(predictions)   
    
def close_session(myrun):
  sess=myrun['sess']
  sess.close()
  return(sess)
    
    
    
    
    

           
    
