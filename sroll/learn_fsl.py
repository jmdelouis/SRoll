import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import sys
#import tensorflow as tf
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()
import time
import netCDF4 as nc4

bolo=['857-1','857-2','857-3','857-4']

dosave=False

if dosave==True:

    rgsize=27664
    nring=26051-240
    for ibolo in bolo:
        print(ibolo)
        ## MODEL
        vv=np.fromfile('/export/home1/jmdeloui/DATA4SROLL4/%s_REP7_2'%(ibolo),offset=240*rgsize*4,count=nring*rgsize,dtype='float32')
        vh=np.fromfile('/export/home1/jmdeloui/DATA4SROLL4/%s_REP6_hit'%(ibolo),offset=240*rgsize*4,count=nring*rgsize,dtype='float32')
        ph=np.fromfile('/export/home1/jmdeloui/DATA4SROLL4/%s_REP6_phregul'%(ibolo),offset=240*rgsize*4,count=nring*rgsize,dtype='float32')
        idx=np.fromfile('/export/home1/jmdeloui/DATA4SROLL4/%s_REP6_rgadutot.int32.bin'%(ibolo),dtype='int32')
        ridx=idx[(np.arange(nring*rgsize)//rgsize)+240]*256+((256*ph/(2*np.pi)).astype('int'))%256

        hh=np.bincount(ridx,weights=vh,minlength=256*256).reshape(256,256)
        rh=np.bincount(ridx,weights=vh*vv,minlength=256*256).reshape(256,256)

        im=rh/hh
        im[hh==0]=np.median(im[hh!=0])

        print(vv[vh>0].std(),(vv-(im.reshape(256*256))[ridx])[vh>0].std())

        #plt.imshow(im,cmap='jet')
        #plt.show()
        np.save('FSL2D_%s_I.npy'%(ibolo),im)
        np.save('FSL2D_%s_H.npy'%(ibolo),hh)
    exit(0)



LEARNING_RATE = 0.3
BATCH_SIZE = 4
DECAY_RATE = 0.9995
NUM_EPOCHS = 4000
EVAL_FREQUENCY = 100

INPARAM = 4

XIMAGE_SIZE=256
YIMAGE_SIZE=256
NCOMP=4
SCALE=2
NDCONV=2
#NDCONV=4
DEEPNESS=5
KERNELSZ=5
NHIDDEN = NCOMP
DOFULLYCONNECTED=True
DORELU=True
XNHIDDEN=XIMAGE_SIZE // (SCALE**NDCONV)
YNHIDDEN=YIMAGE_SIZE // (SCALE**NDCONV)

doloss2=True
NTOTHIDDEN=NHIDDEN*(XNHIDDEN*YNHIDDEN)
NUM_DCONV_CHAN={}

for i in range(NDCONV):
    NUM_DCONV_CHAN[i]=int(DEEPNESS)
NUM_DCONV_CHAN[0] =NHIDDEN
print(XNHIDDEN,YNHIDDEN,NTOTHIDDEN)
NUM_DCONV_CHAN[NDCONV]=1

nfile=1

print("Num GPUs Available: ", len(tf.config.experimental.list_physical_devices('GPU')))
tf.config.set_soft_device_placement(True)
tf.debugging.set_log_device_placement(False)

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

#======================================================================================================================================
#
#======================================================================================================================================
# nout nside de la carte de sortie
# icoef est le poids donné au PWST (10^(icoef)) il semble que 8 soit pas mal. Si on met 0 on suprime l'effet des PWST
# =====================================================================================================================================
# BIDX index du bolometre, chaque bolometre a donc une erreur de calibration
#=====================================================================================================================================

path='/export/home/jmdeloui/857POL/'

def lrelu(x, alpha=0.3):
    #return tf.nn.relu(x)
    return tf.maximum(x, tf.multiply(x, alpha))


res=np.zeros([4,256,256])
resh=np.zeros([4,256,256])
for i in range(4):
    res[i,:,:]=np.load('FSL2D_%s_I.npy'%(bolo[i]))
    resh[i,:,:]=np.load('FSL2D_%s_H.npy'%(bolo[i]))
    
#exit(0)

corrnorm=1/res.std()
data=tf.constant((corrnorm*res).astype('float32'))


def wdecod(parameter):
    # convert into image/bias
    print(fc2_weights.get_shape().as_list())
    print(fc2_biases.get_shape().as_list())
    sys.stderr.write('NPAR PER IMAGE %d\n'%(int((parameter.get_shape().as_list())[1])))
    hidden = tf.matmul(tf.reshape(parameter,[4,INPARAM]), fc2_weights) + fc2_biases
    if DORELU==True:
        hidden = lrelu(hidden)
        
    NTIME=int((parameter.get_shape().as_list())[0])
    sys.stderr.write('NPAR SIZE %d %d\n'%(int((parameter.get_shape().as_list())[0]),int((parameter.get_shape().as_list())[1])))
    
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
                #relu = lrelu(tmp,alpha=0)
                relu=tmp
        else:
            relu=tmp

    relu=relu/100
    return relu


def model(parameter):
    
    modpar=wdecod(parameter)
    return modpar


with tf.device('/gpu:2'):
    fc2_weights = tf.Variable(0.1*np.random.randn(INPARAM,XIMAGE_SIZE // (SCALE**(NDCONV))*YIMAGE_SIZE // (SCALE**(NDCONV))*NUM_DCONV_CHAN[0]).astype('float32'))
    fc2_biases  = tf.Variable(0.1*np.random.randn(XIMAGE_SIZE // (SCALE**(NDCONV))*YIMAGE_SIZE // (SCALE**(NDCONV))*NUM_DCONV_CHAN[0]).astype('float32'))
    dconv_weights={}
    dconv_biases={}
    for i in range(NDCONV):
        dconv_weights[i] = tf.Variable(0.1*np.random.randn(KERNELSZ, KERNELSZ, NUM_DCONV_CHAN[i+1], NUM_DCONV_CHAN[i]).astype('float32'))
        dconv_biases[i] = tf.Variable(0.1*np.random.randn(NUM_DCONV_CHAN[i+1]).astype('float32'))
    
    params = tf.Variable(np.zeros([4,INPARAM],dtype='float32'))

    ampfsl=tf.constant(np.array([[0.35,0.02,0.28,0.15]],dtype='float32'))
    cstfsl=model(params)

    #corrfsl=tf.transpose(model(params))
    corrfsl=tf.reshape(cstfsl,[4,256,256])

with tf.device('/gpu:1'):
    loss=tf.reduce_sum(tf.square(data-corrfsl)) #+1E2*tf.square(tf.reduce_mean(ampfsl)-1.0)

    numbatch = tf.Variable(0, dtype=tf.float32)

    # Decay once per epoch, using an exponential schedule starting at 0.01.
    learning_rate = tf.train.exponential_decay(
        LEARNING_RATE,       # Base learning rate.
        numbatch * BATCH_SIZE,  # Current index into the dataset.
        BATCH_SIZE,          # Decay step.
        DECAY_RATE,          # Decay rate.
        staircase=True)

    # Use simple momentum for the optimization.
    opti=tf.train.AdamOptimizer(learning_rate,0.9).minimize(loss,global_step=numbatch)

# Create a local session to run the training.
start_time = time.time()
start_time_total = start_time

sess = tf.Session(config=tf.ConfigProto(log_device_placement=False,allow_soft_placement=True))
    
# Run all the initializers to prepare the trainable parameters.
tf.global_variables_initializer().run(session=sess)
print('Initialized!')

for step in range(NUM_EPOCHS):
    sess.run(opti)
    if step%EVAL_FREQUENCY==0:
        l, lr,afsl = sess.run([loss, learning_rate,ampfsl])
            
        elapsed_time = time.time() - start_time
        start_time = time.time()
            
        print('Step %d, Mloss[par of sigma]: %7.4f, lrate: %8.5f DURATION %4.2f' % (step,np.sqrt(l/(256*256*4)), lr,elapsed_time),afsl)

idata,thecorr= sess.run([data,corrfsl])

idata/=corrnorm
thecorr/=corrnorm
print(corrnorm)


# -------------------------------------------------------------
#SAVE IN NETCDF

f = nc4.Dataset('FSL_CNN_netcdf.nc','w', format='NETCDF4') 

fsl_grp = f.createGroup('FSL')

#Create dimensions FSL - remove unused dimensions

fsl_grp.createDimension("batchsz",len(bolo))
fsl_grp.createDimension("ximage_size",XIMAGE_SIZE)
fsl_grp.createDimension("yimage_size",YIMAGE_SIZE)
fsl_grp.createDimension("ndconv",NDCONV)
fsl_grp.createDimension("dofullyconnected",float(DOFULLYCONNECTED==True))
fsl_grp.createDimension("deepness",DEEPNESS)
fsl_grp.createDimension("ncomp",NCOMP)
fsl_grp.createDimension("scalex",SCALE)
fsl_grp.createDimension("xkernelsz",KERNELSZ)
fsl_grp.createDimension("scaley",SCALE)
fsl_grp.createDimension("ykernelsz",KERNELSZ)
fsl_grp.createDimension("corrnorm",corrnorm)


## DCONV weights and biases -----
for i in range(NDCONV):
    w,b=sess.run([dconv_weights[i],dconv_biases[i]])
    #create dimensions - A recrire proprement
    fsl_grp.createDimension("f_Xkernel_size_"+str(i),dconv_weights[i].get_shape()[0])
    fsl_grp.createDimension("f_Ykernel_size_"+str(i),dconv_weights[i].get_shape()[1])
    fsl_grp.createDimension("f_Ndin_"+str(i),dconv_weights[i].get_shape()[2])
    fsl_grp.createDimension("f_Ndout_"+str(i),dconv_weights[i].get_shape()[3])
    #Create variables
    fsl_grp.createVariable('dconv_weights_'+str(i), 'f4', ('f_Xkernel_size_'+str(i),'f_Ykernel_size_'+str(i),'f_Ndin_'+str(i),'f_Ndout_'+str(i)))
    fsl_grp.createVariable('dconv_biases_'+str(i), 'f4',('f_Ndin_'+str(i),))
    #load data   
    fsl_grp.variables["dconv_weights_"+str(i)][:,:,:,:]=w
    fsl_grp.variables["dconv_biases_"+str(i)][:,]=b


##  FC weight and biases -----
fw,fb,afsl,apar=sess.run([fc2_weights,fc2_biases,ampfsl,params])
#Create dimensions
fsl_grp.createDimension("dim_mapc2w_0",fw.shape[0])
fsl_grp.createDimension("dim_mapc2w_1",fw.shape[1])
fsl_grp.createDimension("dim_mapc2b",fb.shape[0])
fsl_grp.createDimension("dim_mappar_0",apar.shape[0])
fsl_grp.createDimension("dim_mappar_1",apar.shape[1])
#create variables
fsl_grp.createVariable('mappar', 'f4',('dim_mappar_0','dim_mappar_1'))
fsl_grp.createVariable('fc2_weights', 'f4', ('dim_mapc2w_0','dim_mapc2w_1'))
fsl_grp.createVariable('fc2_biases', 'f4', ('dim_mapc2b',))
#load data
fsl_grp.variables["fc2_weights"][:,:]= fw
fsl_grp.variables["fc2_biases"][:,]= fb
fsl_grp.variables["mappar"][:,:]=apar


#------------------------------------------------------------------------

print('TOTAL DURATION ',time.time() - start_time_total)

for i in range(4):
    
    plt.subplot(3,4,1+i)
    plt.imshow(idata[i,:,:],cmap='jet')
    plt.subplot(3,4,5+i)
    plt.imshow(thecorr[i,:,:],cmap='jet')
    plt.subplot(3,4,9+i)
    plt.imshow(thecorr[i,:,:]-idata[i,:,:],cmap='jet')
plt.show()



