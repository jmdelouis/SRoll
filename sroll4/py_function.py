import numpy as np
import sys
import matplotlib.pyplot as plt
import healpy as hp
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()
import time

TESTCONVERGENCE = False

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

def wdecod(parameter,netinfo,rank,ib):
  # convert into image/bias
  #hidden = tf.matmul(parameter, fc1_weights) + fc1_biases
  #if DORELU==True:
  #    hidden = lrelu(hidden)
  if rank==0:
    sys.stderr.write('NPAR PER IMAGE %d\n'%(int((parameter.get_shape().as_list())[1])))

  if netinfo['DOFULLYCONNECTED']==False:
    hidden = parameter
  else:
    hidden = tf.matmul(parameter, netinfo['fc2_weights'][ib]) + netinfo['fc2_biases'][ib]
    if netinfo['DORELU']==True:
      hidden = lrelu(hidden)
  #hidden = lrelu(hidden)
  NTIME=int((parameter.get_shape().as_list())[0])
  if rank==0:
    sys.stderr.write('NPAR SIZE %d %d\n'%(int((parameter.get_shape().as_list())[0]),int((parameter.get_shape().as_list())[1])))

  relu = tf.reshape(hidden,[NTIME, netinfo['XIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']) , netinfo['YIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']) , netinfo['NUM_DCONV_CHAN'][0]])
  #
  ## deconvolve : First layer
  # DeConvolution kernel
  for i in range(netinfo['NDCONV']):
    conv = tf.nn.conv2d_transpose(relu,netinfo['dconv_weights'][ib][i],strides=[1, netinfo['SCALE'], netinfo['SCALE'], 1],padding='SAME',
                                  output_shape=[NTIME,netinfo['XIMAGE_SIZE'] // (netinfo['SCALE']**(netinfo['NDCONV']-1-i)),netinfo['YIMAGE_SIZE'] // (netinfo['SCALE']**(netinfo['NDCONV']-1-i)),
                                                netinfo['NUM_DCONV_CHAN'][i+1]],name='dconv_%d'%(i))
    tmp=tf.nn.bias_add(conv, netinfo['dconv_biases'][ib][i])
    if rank==0:
      alist=tmp.get_shape().as_list()
      s="["
      for iv in alist:
        s=s+"%d,"%(iv)
      s=s+"],["
      alist=netinfo['dconv_weights'][ib][i].get_shape().as_list()
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
    modpar=wdecod(tf.gather(parameter,0,axis=2),netinfo,rank,0)
    for ib in range(1,netinfo['nbolo']):
        tmp= wdecod(tf.gather(parameter,ib,axis=2),netinfo,rank,ib)
        modpar=tf.concat([modpar,tmp],3)
    return modpar

def alloc_table_float32(ndata):
    mydata = np.zeros([ndata],dtype='float32')
    return mydata

def alloc_table_int32(ndata):
    mydata = np.zeros([ndata],dtype='int32')
    return mydata

def free_table(data):
    del(data)

def init_shape(nbolo,xsize,ysize,rank):
    
    mynetwork={}
    mynetwork['XIMAGE_SIZE']=xsize
    mynetwork['YIMAGE_SIZE']=ysize
    mynetwork['nbolo']=nbolo
    mynetwork['rank']=rank
    mynetwork['NCOMP']=6
    mynetwork['SCALE']=2
    mynetwork['NDCONV']=4
    XNHIDDEN=mynetwork['XIMAGE_SIZE'] // (mynetwork['SCALE']**mynetwork['NDCONV'])
    YNHIDDEN=mynetwork['YIMAGE_SIZE'] // (mynetwork['SCALE']**mynetwork['NDCONV'])
    NHIDDEN=mynetwork['NCOMP']
    DEEPNESS=mynetwork['NCOMP']
    EVAL_FREQUENCY = 100
    mynetwork['KERNELSZ'] = 3

    mynetwork['FromFile']=False
  
    mynetwork['NUM_EPOCHS']= 2000
    SEED = 1234
    mynetwork['DORELU'] = True
    mynetwork['LEARNING_RATE']= 0.03
    mynetwork['DECAY_RATE']=0.995
    mynetwork['DOFULLYCONNECTED'] = False
    mynetwork['BATCHSZ'] = 1
    
    NUM_DCONV_CHAN = {}
    if mynetwork['DOFULLYCONNECTED']==False:
        NTOTHIDDEN=NHIDDEN*(XNHIDDEN*YNHIDDEN)
        for i in range(mynetwork['NDCONV']):
            NUM_DCONV_CHAN[i]=int(DEEPNESS)
        NUM_DCONV_CHAN[0] =NHIDDEN

    else:
        NTOTHIDDEN=1
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

def alloc_param(mynetwork):
  param=np.zeros([mynetwork['BATCHSZ']*mynetwork['NTOTHIDDEN']*mynetwork['nbolo']],dtype='float32')
  return(param)

def get_param(myrun):
  sess=myrun['sess']
  params=sess.run([myrun['params']])[0]
  return(params)

def alloc_convw(mynetwork):
  outobj={}
  for ib in range(mynetwork['nbolo']):
    for i in range(mynetwork['NDCONV']):
      outobj["%03d"%(i+ib*mynetwork['NDCONV'])] = np.zeros([mynetwork['KERNELSZ'] *mynetwork['KERNELSZ'] *mynetwork['NUM_DCONV_CHAN'][i+1]*mynetwork['NUM_DCONV_CHAN'][i]],dtype='float32')
  return(outobj)

def alloc_convb(mynetwork):
  outobj={}
  for ib in range(mynetwork['nbolo']):
    for i in range(mynetwork['NDCONV']):
      outobj["%03d"%(i+ib*mynetwork['NDCONV'])] = np.zeros([mynetwork['NUM_DCONV_CHAN'][i+1]],dtype='float32')
  return(outobj)

def get_convw(myrun):
  sess=myrun['sess']
  outobj={}
  for ib in range(mynetwork['nbolo']):
    for i in range(mynetwork['NDCONV']):
      outobj["%03d"%(i+ib*mynetwork['NDCONV'])] = sess.run([myrun['dconvw']])[0]
  return(outobj)

def get_convb(myrun):
  sess=myrun['sess']
  outobj={}
  for ib in range(mynetwork['nbolo']):
    for i in range(mynetwork['NDCONV']):
      outobj["%03d"%(i+ib*mynetwork['NDCONV'])] = sess.run([myrun['dconvb']])[0]
  return(outobj)

def alloc_flw(netinfo):
  flw={}
  for ib in range(netinfo['nbolo']):
    flw['%03d'%(ib)]=np.zeros([NTOTHIDDEN* netinfo['XIMAGE_SIZE']* (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['YIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['NUM_DCONV_CHAN'][0]],dtype='float32')
  return(flw)

def alloc_flb(mynetwork):
  flb={}
  for ib in range(netinfo['nbolo']):
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
  
def init_network(mynetwork,signal,weights,TCO1,TSI1,MAT0,MAT1,MAT2,hidx,idx,realpix,in_param,in_convw,in_convb,in_flw,in_flb):
    netinfo={}
    for i in mynetwork:
      netinfo[i]=mynetwork[i]
    NUM_DCONV_CHAN=mynetwork['NUM_DCONV_CHAN']
    
    if mynetwork['rank']==0:
      sys.stderr.write('VARIANCE CNN INPUT %f %f\n'%((signal*mynetwork['RAPDATA']).std(),(weights/mynetwork['RAPDATA']).std()))

    NPIX=realpix.shape[0]
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

    coefregul=weights.sum()/mynetwork['RAPDATA']

    params = tf.Variable(in_param.reshape(mynetwork['BATCHSZ'],mynetwork['NTOTHIDDEN'],mynetwork['nbolo']))

    #print(mynetwork['RAPDATA'],signal*mynetwork['RAPDATA'],in_param)
    data   = tf.constant((signal*mynetwork['RAPDATA']).astype('float32'))
    del(signal)
    hdata  = tf.constant((weights/mynetwork['RAPDATA']).astype('float32'))
    del(weights)
    tf_CO  = tf.constant(TCO1)
    del(TCO1)
    tf_SI  = tf.constant(TSI1)
    del(TSI1)
    tf_M0  = tf.constant(MAT0)
    del(MAT0)
    tf_M1  = tf.constant(MAT1)
    del(MAT1)
    tf_M2  = tf.constant(MAT2)
    del(MAT2)
    tf_hidx = tf.constant(hidx)
    del(hidx)
    tf_idx = tf.constant(idx)
    del(idx)


    if netinfo['DOFULLYCONNECTED']:
      if mynetwork['FromFile']==True:
        netinfo['fc2_weights'] = {}
        netinfo['fc2_biases'] = {}
        for ib in range(mynetwork['nbolo']):
          netinfo['fc2_weights'][ib] = tf.constant(mynetwork['fc2w_%d'%(ib)])
          netinfo['fc2_biases'][ib] = tf.constant(mynetwork['fc2b_%d'%(ib)])
      else:
        for ib in range(mynetwork['nbolo']):
          netinfo['fc2_weights_%d'%(ib)] = tf.Variable(in_flw['%03d'%(ib)].reshape(NTOTHIDDEN, netinfo['XIMAGE_SIZE']// (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['YIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['NUM_DCONV_CHAN'][0]))
          netinfo['fc2_biases_%d'%(ib)] = tf.Variable(in_flb['%03d'%(ib)].reshape(netinfo['XIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['YIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['NUM_DCONV_CHAN'][0]))

    dconv_weights={}
    dconv_biases={}
    nw=0
    for ib in range(mynetwork['nbolo']):
      dconv_weights[ib]={}
      dconv_biases[ib]={}
      if mynetwork['FromFile']==True:
        for i in range(netinfo['NDCONV']):
          dconv_weights[ib][i] = tf.constant(mynetwork['lw%d_%d'%(i,ib)])
          dconv_biases[ib][i] = tf.constant(mynetwork["lb%d_%d"%(i,ib)])
      else:
        dconv_weights[ib][i] = tf.Variable(in_convw["%03d_%d"%(i,ib)].reshape(mynetwork['KERNELSZ'],mynetwork['KERNELSZ'],NUM_DCONV_CHAN[i+1],NUM_DCONV_CHAN[i]))
        dconv_biases[ib][i] = tf.Variable(in_convb["%03d_%d"%(i,ib)])

    del(realpix)

    netinfo['dconv_weights']=dconv_weights
    netinfo['dconv_biases']=dconv_biases

    lpar = params
    logits = tf.reshape(model(lpar,netinfo,mynetwork['rank']),[netinfo['BATCHSZ']*netinfo['XIMAGE_SIZE']*netinfo['YIMAGE_SIZE']*netinfo['nbolo']])

    #avvloss=tf.square(tf.reduce_mean(tf.reshape(logits,[netinfo['BATCHSZ'],netinfo['XIMAGE_SIZE']*netinfo['YIMAGE_SIZE']*netinfo['nbolo']]),1))

    vsignal = tf.gather(logits,tf_idx)

    imii =tf.math.unsorted_segment_sum(hdata*vsignal, tf_hidx, NPIX)
    imqq =tf.math.unsorted_segment_sum(hdata*vsignal*tf_CO, tf_hidx, NPIX)
    imuu =tf.math.unsorted_segment_sum(hdata*vsignal*tf_SI, tf_hidx, NPIX)

    residu_vdata = vsignal #- tf_M0*tf.gather(imii,tf_hidx)-tf_M1*tf.gather(imqq,tf_hidx)-tf_M2*tf.gather(imuu,tf_hidx)

    avvloss=tf.square(tf.reduce_sum(hdata*residu_vdata))

    imii =tf.math.unsorted_segment_sum(hdata*data, tf_hidx, NPIX)
    imqq =tf.math.unsorted_segment_sum(hdata*data*tf_CO, tf_hidx, NPIX)
    imuu =tf.math.unsorted_segment_sum(hdata*data*tf_SI, tf_hidx, NPIX)

    #hmap = tf.math.unsorted_segment_sum(hdata*data,tf_idx,netinfo['XIMAGE_SIZE']*netinfo['YIMAGE_SIZE']*netinfo['nbolo'])
    #vmap = tf.math.unsorted_segment_sum(hdata,tf_idx,netinfo['XIMAGE_SIZE']*netinfo['YIMAGE_SIZE']*netinfo['nbolo'])

    residu_data = data #- tf_M0*tf.gather(imii,tf_hidx)-tf_M1*tf.gather(imqq,tf_hidx)-tf_M2*tf.gather(imuu,tf_hidx)

    #loss=tf.reduce_sum(hdata*tf.square(residu_vdata-residu_data)) #+tf.reduce_sum(hdata*tf.square(vsignal))
    #loss = tf.reduce_sum(tf.square(hmap-vmap*logits))
    loss=tf.reduce_sum(hdata*residu_vdata*(residu_vdata-2*residu_data)) #+avvloss

    numbatch = tf.Variable(0, dtype=tf.float32)

    # Decay once per epoch, using an exponential schedule starting at 0.01.
    learning_rate = tf.train.exponential_decay(
        netinfo['LEARNING_RATE'],       # Base learning rate.
        numbatch,            # Current index into the dataset.
        10,                  # Decay step.
        netinfo['DECAY_RATE'],          # Decay rate.
        staircase=True)

    # Use simple momentum for the optimization.
    opti=tf.compat.v1.train.AdamOptimizer(learning_rate,0.9)

    optimizer = opti.minimize(loss,global_step=numbatch)

    igrad={}   
    gradient={}
    assign = {}
    gradient[0] = opti.compute_gradients(loss,var_list=[params])[0]
    igrad[0]    = tf.compat.v1.placeholder(tf.float32,shape=(netinfo['BATCHSZ'],netinfo['NTOTHIDDEN'],netinfo['nbolo']))
    assign[0]   = params.assign(igrad[0])

    nnvar=1

    if mynetwork['FromFile']==False:
      if netinfo['DOFULLYCONNECTED']:
        for ib in range(netinfo['nbolo']):
          igrad[1+ib]    = tf.compat.v1.placeholder(tf.float32,shape=(netinfo['NTOTHIDDEN'], netinfo['XIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['YIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['NUM_DCONV_CHAN'][0]))
          igrad[1+netinfo['nbolo']+ib]     = tf.compat.v1.placeholder(tf.float32,shape=(netinfo['XIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['YIMAGE_SIZE'] // (netinfo['SCALE']**netinfo['NDCONV']) * netinfo['NUM_DCONV_CHAN'][0]))
          gradient[1+ib] = opti.compute_gradients(loss,var_list=[netinfo['fc2_weights']])
          gradient[1+netinfo['nbolo']+ib] = opti.compute_gradients(loss,var_list=[netinfo['fc2_biases']])
        
          assign[1+ib] = netinfo['fc2_weights'].assign(igrad[1+ib])
          assign[1+netinfo['nbolo']+ib] = netinfo['fc2_biases'].assign(igrad[1+netinfo['nbolo']+ib])
        nnvar=1+2*netinfo['nbolo']
      else:
        nnvar=1

      for ib in range(netinfo['nbolo']):
        for i in range(netinfo['NDCONV']):
          igrad[nnvar+i*2*netinfo['nbolo']+ib]    = tf.compat.v1.placeholder(tf.float32,shape=(netinfo['KERNELSZ'], netinfo['KERNELSZ'], NUM_DCONV_CHAN[i+1], NUM_DCONV_CHAN[i]))
          igrad[nnvar+(i*2+1)*netinfo['nbolo']+ib]  = tf.compat.v1.placeholder(tf.float32,shape=(NUM_DCONV_CHAN[i+1]))
        
          gradient[nnvar+i*2*netinfo['nbolo']+ib] = opti.compute_gradients(loss,var_list=[dconv_weights[i]])
          gradient[nnvar+(i*2+1)*netinfo['nbolo']+ib] = opti.compute_gradients(loss,var_list=[dconv_biases[i]])

          assign[nnvar+i*2*netinfo['nbolo']+ib]   = dconv_weights[i].assign(igrad[nnvar+i*2*netinfo['nbolo']+ib])
          assign[nnvar+(i*2+1)*netinfo['nbolo']+ib] = dconv_biases[i].assign(igrad[nnvar+(i*2+1)*netinfo['nbolo']+ib])
          nnvar+=2

    tgradient = [(igrad[i],gradient[i][1]) for i in range(nnvar)] 
    apply_grad = opti.apply_gradients(tgradient,global_step=numbatch)

    if mynetwork['rank']==0:
      sys.stderr.write('Initialized Network\n')

    sess=tf.Session()

    tf.global_variables_initializer().run(session=sess)

    if mynetwork['rank']==0:
      sys.stderr.write('Initialized Variables\n')

    myrun={}
    myrun['loss']=loss
    myrun['vloss']=avvloss
    myrun['regul']=coefregul
    myrun['optimizer']=optimizer
    myrun['sess']=sess
    myrun['logits']=logits
    myrun['learning_rate']=learning_rate
    myrun['params']=params
    if netinfo['DOFULLYCONNECTED']:
        myrun['flw']=netinfo['fc2_weights']
        myrun['flb']=netinfo['fc2_biases']
    myrun['ISFL']=netinfo['DOFULLYCONNECTED']
    myrun['dconvw']=dconv_weights
    myrun['dconvb']=dconv_biases
    myrun['gradient']=gradient
    myrun['igrad']=igrad
    myrun['vsignal']=vsignal
    myrun['assign']=assign
    myrun['apply_grad']=apply_grad
    myrun['rank']=mynetwork['rank']
    myrun['RAPDATA']=mynetwork['RAPDATA']
    myrun['FromFile']=mynetwork['FromFile']
    myrun['nbolo']=netinfo['nbolo']
    myrun['ngradient']=nnvar
    return(myrun)

def calc_grad(myrun):
    sess=myrun['sess']
    vgrad={}
    for i in range(myrun['ngradient']):
      vvv=sess.run([myrun['gradient'][i]])[0]
      vgrad["%03d"%(i)] = np.array(vvv[0])

    return(vgrad)

def calc_opti(myrun):
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

def get_loss(myrun):
    sess=myrun['sess']
    #sess=tf.Session()
    l,lr=sess.run([myrun['loss'],myrun['vloss']])
    res=np.zeros([2],dtype='float32')
    res[0]=l/myrun['regul']
    res[1]=lr/myrun['regul']
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

def get_prediction(myrun):
  sess=myrun['sess']
  predictions=(sess.run([myrun['logits']])[0])/myrun['RAPDATA']
  return(predictions.astype('float32'))

def get_correction(myrun):
  sess=myrun['sess']
  predictions=(sess.run([myrun['vsignal']])[0])/myrun['RAPDATA']
  return(predictions.astype('float32'))

def close_session(myrun):
  sess=myrun['sess']
  sess.close()
  return(sess)
  
           
    
