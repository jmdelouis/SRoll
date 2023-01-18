#!/opt/software/occigen/tools/python/2.7.12/intel/17.0/bin/python

from scipy.ndimage import gaussian_filter
import numpy as np
import sys
#import tensorflow as tf
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()
import time
import logging
import os
import argparse
import matplotlib.pyplot as plt

DATAPATH='/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/PBR_JMD'
OUTNAME='545_FSL'
OUTNAME='857_FSL'

boloref=[sys.argv[1]]
TFLEARN=DATAPATH+'/'+'FSL_545_INITMODEL'

dostep=int(sys.argv[2])
ostep=dostep
TFLEARN=DATAPATH+'/'+'FSL_%s_S%d_INITMODEL'%(boloref[0],ostep)
print(TFLEARN)

do857=True
doDef=False
docor=False
dovar=False
doadd=False
doplot=False

if len(sys.argv)>3:
    doplot=True
    
lamptab=[27.5,2.0,20.0,18.0]

if sys.argv[1]=='857-1':
    lamp=27.5
    positivity=1.0
    #======================================================================
    #=    
    #=    REGLAGE THRESHOLD 857-1
    #=    
    #======================================================================
    thresh=0.0007
if sys.argv[1]=='857-2':
    positivity=1.0
    lamp=2.0
    thresh=0.02
if sys.argv[1]=='857-3':
    positivity=-1.0
    lamp=20.0
    thresh=0.0012
if sys.argv[1]=='857-4':
    positivity=1.0
    lamp=8.0
    thresh=0.002

if do857==False:
    res=np.load(DATAPATH+'/INITCODE_%s_ALL.npy'%(OUTNAME))
    nx,ny,nz=res.shape
    res=np.append(np.load(DATAPATH+'/INITCODE_%s_ALL.npy'%(OUTNAME)),np.load(DATAPATH+'/INITCODE1_%s_ALL.npy'%(OUTNAME))).reshape(2*nx,ny,nz)
    hres=np.append(np.load(DATAPATH+'/INITCODE_%s_HALL.npy'%(OUTNAME)),np.load(DATAPATH+'/INITCODE1_%s_HALL.npy'%(OUTNAME))).reshape(2*nx,ny,nz)
else:
    #res=np.load(DATAPATH+'/INITCODE_%s_ALL.npy'%(OUTNAME))[0:4,:,:]
    #hres=np.load(DATAPATH+'/INITCODE_%s_HALL.npy'%(OUTNAME))[0:4,:,:]
    if dostep>0:
        ires=np.load(DATAPATH+'/ERRCODE_%s_%s_ALL.npy'%(OUTNAME,sys.argv[1]))
        ihres=np.load(DATAPATH+'/ERRCODE_%s_%s_HALL.npy'%(OUTNAME,sys.argv[1]))
    
        nx,ny,nz=ires.shape
        bolo=['857-1','857-2','857-3','857-4']
        tmp=np.load('TMPFSL.npy')
        res=np.zeros([5,ny,nz])
        hres=np.zeros([5,ny,nz])
        res[0,:,:]=ires[0,:,:]
        hres[0,:,:]=ihres[0,:,:]
        if doplot==True:
            tmp=np.load('TMPFSL.npy')
            plt.subplot(1,5,1)
            plt.imshow(tmp[0,:,:]+tmp[1,:,:]+tmp[2,:,:]+tmp[3,:,:],cmap='jet')
        
        for i in range(4):
            res2=np.load(DATAPATH+'/ERRCODE_%s_%s_ALL.npy'%(OUTNAME,bolo[i]))
            hres2=np.load(DATAPATH+'/ERRCODE_%s_%s_HALL.npy'%(OUTNAME,bolo[i]))
            if dostep==1:
                mm=np.load('/export/home1/jmdeloui/DATA4SROLL4/VECT/%s_offsets_PROD_REP6_RD12_545GHz_DATA_CNN2D_TL_FSL_1_-32_DECOD_0.npy'%(bolo[i])).reshape(1024,512)
            else:
                mm=np.load('/export/home1/jmdeloui/DATA4SROLL4/VECT/%s_offsets_PROD_REP6_RD12_545GHz_DATA_CNN2D_TL_FSL_S%d_1_-32_DECOD_0.npy'%(bolo[i],dostep-1)).reshape(1024,512)

            aa=((mm*res2[0,:,:])*(hres2[0,:,:]>0)*(hres[0,:,:]>0)).sum()/((mm*mm)*(hres2[0,:,:]>0)*(hres[0,:,:]>0)).sum()
            print(aa)
            aa=1E-3
            template=2*((i!=2)-0.5)*(mm*aa-res2[0,:,:])
            template*=(abs(template)>0.01) #>0.01*abs(template).max())
            res[i+1,:,:]=(ires[0,:,:]+positivity*template)*(hres2[0,:,:]>0)*(hres[0,:,:]>0)
            hres[i+1,:,:]=ihres[0,:,:]*(hres2[0,:,:]>0)*(hres[0,:,:]>0)
            if doplot==True:
                plt.subplot(1,5,2+i)
                plt.imshow((res[i+1,:,:]-ires[0,:,:])*(hres2[0,:,:]>0)*(hres[0,:,:]>0),cmap='jet')
        if doplot==True:
            plt.show()
            exit(0)

    else:
        #======================================================================
        #=    ON NE PASSE QU'ICI
        #=    ON NE PASSE QU'ICI
        #=    ON NE PASSE QU'ICI
        #======================================================================
        res=np.load(DATAPATH+'/ERRCODE_%s_%s_ALL.npy'%(OUTNAME,sys.argv[1]))
        hres=np.load(DATAPATH+'/ERRCODE_%s_%s_HALL.npy'%(OUTNAME,sys.argv[1]))
        nx,ny,nz=res.shape
            
        if dostep<0:
            res=np.load(DATAPATH+'/ERRCODE_%s_%s_ALL.npy'%(OUTNAME,'857-2'))
            hres=np.load(DATAPATH+'/ERRCODE_%s_%s_HALL.npy'%(OUTNAME,'857-2'))
            nx,ny,nz=res.shape

            tmp=np.zeros([3,ny,nz])
            htmp=np.zeros([3,ny,nz])
            thresh=0.001
            for i in range(3):
                if i>=1:
                    ii=i+1
                else:
                    ii=i
                print(lamptab[ii])
                tmp[i,:,:]=(res[1+i,:,:])*(abs(res[0,:,:])>thresh)/lamptab[ii]
                htmp[i,:,:]=hres[1+i,:,:]
                plt.subplot(1,3,1+i)
                #=================================================================================================
                # SUPPRIME LES STRUCTURES NON LIEES AU FSL DANS LE MODEL 3-2
                if i==1:
                    vv=1E-3
                    oo=4E-5
                    ww=(np.repeat(np.exp(-vv*(np.arange(ny)-480)**2),nz)*np.tile(np.exp(-oo*(np.arange(nz)-50)**2),ny)).reshape(ny,nz)
                    ww+=(np.repeat(np.exp(-vv*(np.arange(ny)-30)**2),nz)*np.tile(np.exp(-oo*(np.arange(nz)-50)**2),ny)).reshape(ny,nz) 
                    ww+=(np.repeat(np.exp(-vv*(np.arange(ny)-930)**2),nz)*np.tile(np.exp(-oo*(np.arange(nz)-50)**2),ny)).reshape(ny,nz) 
                    ww+=(np.repeat(np.exp(-vv*(np.arange(ny)-670)**2),nz)*np.tile(np.exp(-oo*(np.arange(nz)-460)**2),ny)).reshape(ny,nz) 
                    ww+=(np.repeat(np.exp(-vv*(np.arange(ny)-220)**2),nz)*np.tile(np.exp(-oo*(np.arange(nz)-460)**2),ny)).reshape(ny,nz)
                    tmp[i,:,:]=(1-ww)*tmp[1,:,:]
                
                plt.imshow(tmp[i,:,:],cmap='jet')
                plt.ylabel('res%d'%(i))
            plt.show()

            print('save the result')
            np.save('TMPFSL2.npy',tmp)
            np.save('TMPFSL2H.npy',htmp)
            exit(0)
        
        tmp=np.load('TMPFSL2.npy')
        htmp=np.load('TMPFSL2H.npy')
        for i in range(3):
            res[1+i,:,:]=res[0,:,:]+tmp[i,:,:]
            hres[1+i,:,:]=hres[0,:,:]*(htmp[i,:,:]>0)
        dobolo2=False
        if dobolo2==True:
            #============================================================================
            # ONLY LEARN DIFF WITH BOLOMETER 2
            #============================================================================
            if sys.argv[1]=='857-1':
                for k in range(3):
                    res[1+k:,:]=res[0,:,:]+tmp[0,:,:]*(3*k-1.0)
                    hres[1+k,:,:]=hres[0,:,:]*(htmp[0,:,:]>0)
            if sys.argv[1]=='857-2':
                nx,ny,nz=res.shape
                for k in range(3):
                    res[1+k,:,:]=res[0,:,:]+np.random.randn(ny,nz)*0.1*res[0,:,:].std()
            if sys.argv[1]=='857-3':
                for k in range(3):
                    res[1+k,:,:]=res[0,:,:]+tmp[1,:,:]*(3*k-1.0)
                    hres[1+k,:,:]=hres[0,:,:]*(htmp[1,:,:]>0)
            if sys.argv[1]=='857-4':
                for k in range(3):
                    res[1+k,:,:]=res[0,:,:]-tmp[2,:,:]*(3*k-1.0)
                    hres[1+k,:,:]=hres[0,:,:]*(htmp[2,:,:]>0)
        
        if doplot==True:
            plt.subplot(1,4,1)
            plt.imshow(res[0,:,:],cmap='jet',vmin=-0.1*res[0,:,:].max(),vmax=res[0,:,:].max())
            for i in range(3):
                plt.subplot(1,4,2+i)
                plt.imshow(res[1+i,:,:]-res[0,:,:],cmap='jet',vmin=-0.1*res[0,:,:].max(),vmax=0.1*res[0,:,:].max())
            plt.show()
            exit(0)
            
    res=res.astype('float32')
    hres=hres.astype('float32')
        

nx,ny,nz=res.shape

if doadd==True:
    temp=np.load(DATAPATH+'/ERRCODE_%s_ALL.npy'%(OUTNAME))
    htemp=np.load(DATAPATH+'/ERRCODE_%s_HALL.npy'%(OUTNAME))
    ncorr,ny,nz=temp.shape
    tmp  = np.zeros([(nx+ncorr*2),ny,nz])
    htmp = np.zeros([(nx+ncorr*2),ny,nz])

    for i in range(nx):
        tmp[i,:,:]=res[i,:,:]
        htmp[i,:,:]=hres[i,:,:]

    ttt=[0,2,3]
    
    for j in range(ncorr):
        tmp[nx+2*j,:,:]    = res[ttt[j],:,:]+temp[j,:,:]/8.0
        htmp[nx+2*j,:,:]   = hres[ttt[j],:,:]
        tmp[nx+2*j+1,:,:]  = res[ttt[j],:,:]-temp[j,:,:]/8.0
        htmp[nx+2*j+1,:,:] = hres[ttt[j],:,:]
        
    res=tmp.astype('float32')
    hres=htmp.astype('float32')
    
print(res.shape)
nx,ny,nz=res.shape
for i in range(nx-1):
    res[i+1,:,:]=res[i+1,:,:]/res[i+1,:,:].std()*res[i,:,:].std()


nd,ny,nx=res.shape

if dovar==True:
    nharm=1
    tmp=np.zeros([(1+4*nharm)*nd,ny,nx])
    htmp=np.zeros([(1+4*nharm)*nd,ny,nx])
    amp=0.2
    
    x = np.arange(nx)
    y = np.arange(ny)
    xx, yy = np.meshgrid(x, y)
    
    for i in range(nd):
        
        tmp[i*(1+4*nharm),:,:]=res[i,:,:]
        htmp[i*(1+4*nharm),:,:]=hres[i,:,:]
        plt.subplot(nd,1+4*nharm,1+i*(1+4*nharm))
        plt.imshow(tmp[i*(1+4*nharm),:,:],cmap='jet')
        for j in range(nharm):
            tmp[i*(1+4*nharm)+1+4*j,:,:]=res[i,:,:]*(1+amp*np.sin(xx/nx*2*np.pi*(1+j)))
            htmp[i*(1+4*nharm)+1+4*j,:,:]=hres[i,:,:]
            plt.subplot(nd,1+4*nharm,2+4*j+i*(1+4*nharm))
            plt.imshow(tmp[i*(1+4*nharm)+1+4*j,:,:],cmap='jet')
            tmp[i*(1+4*nharm)+1+4*j+1,:,:]=res[i,:,:]*(1-amp*np.sin(xx/nx*2*np.pi*(1+j)))
            htmp[i*(1+4*nharm)+1+4*j+1,:,:]=hres[i,:,:]
            plt.subplot(nd,1+4*nharm,3+4*j+i*(1+4*nharm))
            plt.imshow(tmp[i*(1+4*nharm)+1+4*j+1,:,:],cmap='jet')
            tmp[i*(1+4*nharm)+1+4*j+2,:,:]=res[i,:,:]*(1+amp*np.sin(yy/ny*2*np.pi*(1+j)))
            htmp[i*(1+4*nharm)+1+4*j+2,:,:]=hres[i,:,:]
            plt.subplot(nd,1+4*nharm,4+4*j+i*(1+4*nharm))
            plt.imshow(tmp[i*(1+4*nharm)+1+4*j+2,:,:],cmap='jet')
            tmp[i*(1+4*nharm)+1+4*j+3,:,:]=res[i,:,:]*(1-amp*np.sin(yy/ny*2*np.pi*(1+j)))
            htmp[i*(1+4*nharm)+1+4*j+3,:,:]=hres[i,:,:]
            plt.subplot(nd,1+4*nharm,5+4*j+i*(1+4*nharm))
            plt.imshow(tmp[i*(1+4*nharm)+1+4*j+3,:,:],cmap='jet')
    plt.show()
    res=tmp.astype('float32')
    hres=htmp.astype('float32')
    
nd,ny,nx=res.shape

if docor==True:
    ntest=4
    
    x0=341
    y0a=343
    y0b=793

    x1=90
    y1a=386
    y1b=835
    
    sx=0.0006
    sy=0.03

    tmp=np.zeros([(1+2*ntest)*nd,ny,nx])
    htmp=np.zeros([(1+2*ntest)*nd,ny,nx])
    
    x = np.arange(nx)
    y = np.arange(ny)
    xx, yy = np.meshgrid(x, y)
    for i in range(nd):
        dx=30*np.random.randn(ntest)
        dy=30*np.random.randn(ntest)
        
        tmp[i*(1+2*ntest),:,:]=res[i,:,:]
        htmp[i*(1+2*ntest),:,:]=hres[i,:,:]
        plt.subplot(nd,1+2*ntest,1+i*(1+2*ntest))
        plt.imshow(tmp[i*(1+2*ntest),:,:],cmap='jet')
        for j in range(ntest):
            tmp[i*(1+2*ntest)+1+2*j,:,:]=res[i,:,:]+0.01*(np.exp(-(sx*(xx-x0+dx[j])**2+sy*(yy-y0a+dy[j])**2))+np.exp(-(sx*(xx-x0+dx[j])**2+sy*(yy-y0b+dy[j])**2)))
            htmp[i*(1+2*ntest)+1+2*j,:,:]=hres[i,:,:]
            plt.subplot(nd,1+2*ntest,2+2*j+i*(1+2*ntest))
            plt.imshow(tmp[i*(1+2*ntest)+1+2*j,:,:],cmap='jet')
            tmp[i*(1+2*ntest)+1+2*j+1,:,:]=res[i,:,:]+0.01*(np.exp(-(sx*(xx-x1+dx[j])**2+sy*(yy-y1a+dy[j])**2))+np.exp(-(sx*(xx-x1+dx[j])**2+sy*(yy-y1b+dy[j])**2)))
            htmp[i*(1+2*ntest)+1+2*j+1,:,:]=hres[i,:,:]
            plt.subplot(nd,1+2*ntest,3+2*j+i*(1+2*ntest))
            plt.imshow(tmp[i*(1+2*ntest)+1+2*j+1,:,:],cmap='jet')
    plt.show()
    res=tmp.astype('float32')
    hres=htmp.astype('float32')
        
import scipy.special as sp

def calcshift(nx,sig=128.0,phase=0.0,amp=0.5):
    xx=np.ones([nx])+np.exp(-0.5*((np.arange(nx)-phase)/sig)**2)
    xx=np.cumsum(xx)
    xx-=xx.min()
    xx=nx*xx/xx.max()
    return(xx)

nd,ny,nx=res.shape
if doDef==True:
    tmp=np.zeros([3*nd,ny,nx])
    htmp=np.zeros([3*nd,ny,nx])

    import scipy.interpolate as interpolate

    x = np.arange(nx)
    y = np.arange(ny)
    xx, yy = np.meshgrid(x, y)
    for i in range(nd):
        f = interpolate.interp2d(x, y, res[i,:,:], kind='cubic')
        fh = interpolate.interp2d(x, y, hres[i,:,:], kind='cubic')
        tmp[i*3,:,:]=res[i,:,:]
        htmp[i*3,:,:]=hres[i,:,:]
        tmp[i*3+1,:,:]=f(calcshift(nx,phase=nx/4),y)
        htmp[i*3+1,:,:]=fh(calcshift(nx,phase=nx/4),y)
        tmp[i*3+2,:,:]=f(calcshift(nx,phase=3*nx/4),y)
        htmp[i*3+2,:,:]=fh(calcshift(nx,phase=3*nx/4),y)
        
    res=tmp.astype('float32')
    hres=htmp.astype('float32')

nx,ny,nz=res.shape
for i in range(nx-1):
    res[i+1,:,:]=res[i+1,:,:]/res[i+1,:,:].max()*res[i,:,:].max()

#================================================================================
#
#    MULTIPLY BY 1000 TO MAKE IT AT THE GOOD LEVEL (x1000)
#
#================================================================================

res=res*1000.
XNHIDDEN=2
YNHIDDEN=XNHIDDEN
NCOMP=16
SCALE=4
NHIDDEN=XNHIDDEN*YNHIDDEN*NCOMP
DEEPNESS=2*NCOMP
EVAL_FREQUENCY = 100
KERNELSZ = 7
SEED = 1234
DORELU = True
LEARNING_RATE= 0.03
DECAY_RATE=0.995
NUM_EPOCHS= 10000
BATCHSZ=nx
DOFULLYCONNECTED = True

if DOFULLYCONNECTED==False:
    NTOTHIDDEN=NHIDDEN//(XNHIDDEN*YNHIDDEN)
    NUM_DCONV_CHAN =[ NHIDDEN//(XNHIDDEN*YNHIDDEN), DEEPNESS, DEEPNESS, 1]
else:
    NTOTHIDDEN=1
    NUM_DCONV_CHAN =[ DEEPNESS, DEEPNESS, DEEPNESS, 1]
NDCONV=len(NUM_DCONV_CHAN)-1
NMAT,XIMAGE_SIZE,YIMAGE_SIZE=res.shape

if DOFULLYCONNECTED:
    fc2_weights = tf.Variable(0.1*np.random.randn(NTOTHIDDEN, XIMAGE_SIZE // (SCALE**NDCONV) * YIMAGE_SIZE // (SCALE**NDCONV) * NUM_DCONV_CHAN[0]).astype('float32'))
    fc2_biases = tf.Variable(0.1*np.random.randn(XIMAGE_SIZE // (SCALE**NDCONV) * YIMAGE_SIZE // (SCALE**NDCONV) * NUM_DCONV_CHAN[0]).astype('float32'))

dconv_weights={}
dconv_biases={}
nw=0
for i in range(NDCONV):
    dconv_weights[i] = tf.Variable(0.1*np.random.randn(KERNELSZ, KERNELSZ, NUM_DCONV_CHAN[i+1], NUM_DCONV_CHAN[i]).astype('float32'))
    dconv_biases[i] = tf.Variable(0.1*np.random.randn(NUM_DCONV_CHAN[i+1]).astype('float32'))

params = tf.Variable(np.zeros([NMAT,NTOTHIDDEN],dtype='float32'))

train_data = tf.placeholder(tf.float32,shape=(BATCHSZ,XIMAGE_SIZE,YIMAGE_SIZE,1))
hdata=tf.constant(hres.reshape(BATCHSZ,XIMAGE_SIZE,YIMAGE_SIZE,1))
run_param  = tf.placeholder(tf.float32,shape=(1,NTOTHIDDEN))


def lrelu(x, alpha=0.3):
  #return tf.nn.relu(x)
  return tf.maximum(x, tf.multiply(x, alpha))
    
def wdecod(parameter):
    # convert into image/bias
    #hidden = tf.matmul(parameter, fc1_weights) + fc1_biases
    #if DORELU==True:
    #    hidden = lrelu(hidden)
    print('NPAR PER IMAGE',int((parameter.get_shape().as_list())[1]))
    if DOFULLYCONNECTED==False:
        hidden = parameter
    else:
        hidden = tf.matmul(parameter, fc2_weights) + fc2_biases
        if DORELU==True:
            hidden = lrelu(hidden)
    #hidden = lrelu(hidden)
    NTIME=int((parameter.get_shape().as_list())[0])

    relu = tf.reshape(hidden,[NTIME, XIMAGE_SIZE // (SCALE**NDCONV) , YIMAGE_SIZE // (SCALE**NDCONV) , NUM_DCONV_CHAN[0]])
    #
    ## deconvolve : First layer
    # DeConvolution kernel
    for i in range(NDCONV):
        conv = tf.nn.conv2d_transpose(relu,dconv_weights[i],strides=[1, SCALE, SCALE, 1],padding='SAME',
                                      output_shape=[NTIME,XIMAGE_SIZE // (SCALE**(NDCONV-1-i)),YIMAGE_SIZE // (SCALE**(NDCONV-1-i)),
                                                    NUM_DCONV_CHAN[i+1]],name='dconv_%d'%(i))
        tmp=tf.nn.bias_add(conv, dconv_biases[i])
        print(tmp.get_shape().as_list(),dconv_weights[i].get_shape().as_list())
        
        if DORELU==True:
            if i!=NDCONV-1:
                relu = lrelu(tmp)
            else:
                relu = lrelu(tmp,alpha=0)
        else:
            relu=tmp

    relu=relu/100
      
    return relu
  

def model(parameter):

    modpar=wdecod(parameter)
    return modpar


lpar = params
logits = model(lpar)

synthe = model(run_param)

loss=tf.reduce_sum(hdata*tf.square((logits-train_data)))


numbatch = tf.Variable(0, dtype=tf.float32)
    
# Decay once per epoch, using an exponential schedule starting at 0.01.
learning_rate = tf.train.exponential_decay(
        LEARNING_RATE,       # Base learning rate.
        numbatch,            # Current index into the dataset.
        10,                  # Decay step.
        DECAY_RATE,          # Decay rate.
        staircase=True)
    
    
optimizer=tf.train.AdamOptimizer(learning_rate,0.9).minimize(loss,
                                                                           global_step=numbatch)
# Create a local session to run the training.
start_time = time.time()
print('Initialized 1')

import matplotlib.pyplot as plt

loss_log=np.zeros([NUM_EPOCHS//EVAL_FREQUENCY])

minloss=1E30
with tf.Session() as sess:
    # Run all the initializers to prepare the trainable parameters.
    tf.global_variables_initializer().run()
    print('Initialized 2')
    # Loop through training steps.
      
    for step in range(NUM_EPOCHS):
        
        feed_dict = {train_data: res.reshape(BATCHSZ,XIMAGE_SIZE,YIMAGE_SIZE,1)}
      
        sess.run(optimizer, feed_dict=feed_dict)
        # print some extra information once reach the evaluation frequency
        if step % EVAL_FREQUENCY == 0:
            # fetch some extra nodes' data
            l, lr = sess.run([loss, learning_rate],feed_dict=feed_dict)
            loss_log[step//EVAL_FREQUENCY]=l
                
            elapsed_time = time.time() - start_time
            start_time = time.time()
        
            print('Step %d[%5.2f], Mloss[par of sigma]: %7.4f, lrate: %8.5f DURATION %4.2f' % (step,step*BATCHSZ/(1.0*NMAT),np.sqrt(l/(BATCHSZ*XIMAGE_SIZE*YIMAGE_SIZE)), lr,elapsed_time))
    
                
    lf,lb =sess.run([fc2_weights,fc2_biases])
    np.save('%s_%d_fc2w.npy'%(TFLEARN,0),lf)
    np.save('%s_%d_fc2b.npy'%(TFLEARN,0),lb)
    
    for i in range(NDCONV):
        lw,lb =sess.run([dconv_weights[i],dconv_biases[i]])
        np.save('%s_%d_w_%d.npy'%(TFLEARN,0,i),lw)
        np.save('%s_%d_b_%d.npy'%(TFLEARN,0,i),lb)

    predictions,par=sess.run([logits,lpar], feed_dict=feed_dict)
    np.save('%s_DECODREBUILD.npy'%(TFLEARN),predictions.reshape(NMAT,XIMAGE_SIZE,YIMAGE_SIZE))
    alldata=predictions
    allparams=par
    for i in range(BATCHSZ):
        np.save('%s_%d_param.npy'%(TFLEARN,i),par[i,:].reshape(1,NTOTHIDDEN))
        minloss=l
        np.save('%s_step.npy'%(TFLEARN),ostep)
        print('MIN SAVE',l)
    print(allparams)
    
