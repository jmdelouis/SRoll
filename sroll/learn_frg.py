import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import sys
#import tensorflow as tf
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()
import time
import netCDF4 as nc4

nout=32
dosave=False
if dosave==True:
   mapq=hp.ud_grade(np.fromfile('/export/home/tfoulquier/data_sroll/MAP/353_Q.float32.bin',dtype='float32'),nout)
   mapu=hp.ud_grade(np.fromfile('/export/home/tfoulquier/data_sroll/MAP/353_U.float32.bin',dtype='float32'),nout)

   idx=hp.nest2ring(nout,np.arange(12*nout*nout))
   hp.mollview(mapq[idx],cmap='jet',norm='hist',nest=True)
   hp.mollview(mapu[idx],cmap='jet',norm='hist',nest=True)
   plt.show()
   res=np.zeros([12*nout*nout,2])
   res[:,0]=mapq[idx]
   res[:,1]=mapu[idx]
   np.save('FRG_DATA.npy',res)
res=np.load('FRG_DATA.npy')

LEARNING_RATE = 0.3
BATCH_SIZE = 4
DECAY_RATE = 0.9995
NUM_EPOCHS = 4000
EVAL_FREQUENCY = 100

BATCHSZ=1
INPARAM = 4
XIMAGE_SIZE=12*nout*nout
YIMAGE_SIZE=1
NCOMP=4
SCALEX=4
SCALEY=1
NDCONV=2
#NDCONV=1
DEEPNESS=12
XKERNELSZ=4
YKERNELSZ=1
NHIDDEN = NCOMP
DOFULLYCONNECTED=True
DORELU=True
XNHIDDEN=XIMAGE_SIZE // (SCALEX**NDCONV)
YNHIDDEN=YIMAGE_SIZE // (SCALEY**NDCONV)
if DOFULLYCONNECTED==False:
    INPARAM=XIMAGE_SIZE // (SCALEX**NDCONV) *YIMAGE_SIZE // (SCALEY**NDCONV)* NCOMP

doloss2=True
NTOTHIDDEN=NHIDDEN*(XNHIDDEN*YNHIDDEN)
NUM_DCONV_CHAN={}

for i in range(NDCONV):
    NUM_DCONV_CHAN[i]=int(DEEPNESS)
NUM_DCONV_CHAN[0] = NHIDDEN
print(XNHIDDEN,YNHIDDEN,NTOTHIDDEN)
NUM_DCONV_CHAN[NDCONV]=2

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


#exit(0)

corrnorm=1/res.std()
data=tf.constant((corrnorm*res).astype('float32'))



def wdecod(parameter):
    # convert into image/bias
    if DOFULLYCONNECTED==True:
        print(fc2_weights.get_shape().as_list())
        print(fc2_biases.get_shape().as_list())
        sys.stderr.write('NPAR PER IMAGE %d\n'%(int((parameter.get_shape().as_list())[1])))
        hidden = tf.matmul(tf.reshape(parameter,[1,INPARAM]), fc2_weights) + fc2_biases
    else:
        hidden = parameter
    if DORELU==True:
        hidden = lrelu(hidden)
        
    NTIME=int((parameter.get_shape().as_list())[0])
    sys.stderr.write('NPAR SIZE %d %d\n'%(int((parameter.get_shape().as_list())[0]),int((parameter.get_shape().as_list())[1])))
    
    relu = tf.reshape(hidden,[NTIME, XIMAGE_SIZE // (SCALEX**NDCONV) , YIMAGE_SIZE // (SCALEY**NDCONV) , NUM_DCONV_CHAN[0]])
    #
    ## deconvolve : First layer
    # DeConvolution kernel
    for i in range(NDCONV):
        conv = tf.nn.conv2d_transpose(relu,dconv_weights[i],strides=[1, SCALEX, SCALEY, 1],padding='SAME',
                                      output_shape=[NTIME,XIMAGE_SIZE // (SCALEX**(NDCONV-1-i)),YIMAGE_SIZE // (SCALEY**(NDCONV-1-i)),
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

if DOFULLYCONNECTED:
    fc2_weights = tf.Variable(0.1*np.random.randn(INPARAM,XIMAGE_SIZE // (SCALEX**(NDCONV))*YIMAGE_SIZE // (SCALEY**(NDCONV))*NUM_DCONV_CHAN[0]).astype('float32'))
    fc2_biases  = tf.Variable(0.1*np.random.randn(XIMAGE_SIZE // (SCALEX**(NDCONV))*YIMAGE_SIZE // (SCALEY**(NDCONV))*NUM_DCONV_CHAN[0]).astype('float32'))
dconv_weights={}
dconv_biases={}
for i in range(NDCONV):
    dconv_weights[i] = tf.Variable(0.1*np.random.randn(XKERNELSZ, YKERNELSZ, NUM_DCONV_CHAN[i+1], NUM_DCONV_CHAN[i]).astype('float32'))
    dconv_biases[i] = tf.Variable(0.1*np.random.randn(NUM_DCONV_CHAN[i+1]).astype('float32'))
    
params = tf.Variable(np.zeros([BATCHSZ,INPARAM],dtype='float32'))

ampfsl=tf.constant(np.array([[0.35,0.02,0.28,0.15]],dtype='float32'))
cstfsl=model(params)

#corrfsl=tf.transpose(model(params))
corrfsl=tf.reshape(cstfsl,[12*nout*nout,2])

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

#======================================================================================
## SAVED IN NETCDF
f = nc4.Dataset('MAP_CNN_netcdf.nc','w', format='NETCDF4') 

map_grp = f.createGroup('MAP')

#Create dimensions MAP
map_grp.createDimension("batchsz",BATCHSZ)
map_grp.createDimension("ximage_size",XIMAGE_SIZE)
map_grp.createDimension("yimage_size",YIMAGE_SIZE)
map_grp.createDimension("ndconv",NDCONV)
map_grp.createDimension("dofullyconnected",float(DOFULLYCONNECTED==True))
map_grp.createDimension("deepness",DEEPNESS)
map_grp.createDimension("ncomp",NCOMP)
map_grp.createDimension("scalex",SCALEX)
map_grp.createDimension("xkernelsz",XKERNELSZ)
map_grp.createDimension("scaley",SCALEY)
map_grp.createDimension("ykernelsz",YKERNELSZ)
map_grp.createDimension("corrnorm",corrnorm)

#FC weight and biases
if DOFULLYCONNECTED==True:
  fw,fb,afsl,apar=sess.run([fc2_weights,fc2_biases,ampfsl,params])
  #Create dimensions
  map_grp.createDimension("dim_mapc2w_0",fw.shape[0])
  map_grp.createDimension("dim_mapc2w_1",fw.shape[1])
  map_grp.createDimension("dim_mapc2b",fb.shape[0])
  map_grp.createDimension("dim_mappar_0",apar.shape[0])
  map_grp.createDimension("dim_mappar_1",apar.shape[1])
  #Create variables
  map_grp.createVariable('fc2_weights', 'f4', ('dim_mapc2w_0','dim_mapc2w_1'))
  map_grp.createVariable('fc2_biases', 'f4', ('dim_mapc2b',))
  map_grp.createVariable('mappar', 'f4',('dim_mappar_0','dim_mappar_1'))
  #load data : weight and biases fullyconnected
  map_grp.variables["fc2_weights"][:,:]=fw
  map_grp.variables["fc2_biases"][:,]=fb
  map_grp.variables["mappar"][:,:]=apar

else:
    afsl,apar=sess.run([ampfsl,params])
    #create dimensions
    map_grp.createDimension("dim_mappar_0",apar.shape[0])
    map_grp.createDimension("dim_mappar_1",apar.shape[1])
    map_grp.createDimension("dim_mapamp_0",apar.shape[0])
    map_grp.createDimension("dim_mapamp_1",apar.shape[1])
    #Create variables
    map_grp.createVariable('mappar', 'f4',('dim_mappar_0','dim_mappar_1'))
    map_grp.createVariable('mapamp', 'f4',('dim_mapamp_0','dim_mapamp_1'))
    #load data
    map_grp.variables["mappar"][:,:]=apar
    map_grp.variables["mapamp"][:,:]=afsl

## DCONV weights and biases
for i in range(NDCONV):
    #create dimensions - A recrire proprement
    map_grp.createDimension("m_Xkernel_size_"+str(i),dconv_weights[i].get_shape()[0])
    map_grp.createDimension("m_Ykernel_size_"+str(i),dconv_weights[i].get_shape()[1])
    map_grp.createDimension("m_Ndin_"+str(i),dconv_weights[i].get_shape()[2])
    map_grp.createDimension("m_Ndout_"+str(i),dconv_weights[i].get_shape()[3])
    #Create variables
    map_grp.createVariable('dconv_weights_'+str(i), 'f4', ('m_Xkernel_size_'+str(i),'m_Ykernel_size_'+str(i),'m_Ndin_'+str(i),'m_Ndout_'+str(i)))
    map_grp.createVariable('dconv_biases_'+str(i), 'f4',('m_Ndin_'+str(i),))
    #load data
    w,b=sess.run([dconv_weights[i],dconv_biases[i]])
    map_grp.variables["dconv_weights_"+str(i)][:,:,:,:]=w
    map_grp.variables["dconv_biases_"+str(i)][:,]=b

f.close()


#======================================================================================
print('TOTAL DURATION ',time.time() - start_time_total)

idata/=corrnorm
thecorr/=corrnorm
print(corrnorm)
np.save('/export/home/tfoulquier/heal_new/MAPAMP.npy',afsl)
np.save('/export/home/tfoulquier/heal_new/MAPDATA.npy',idata)
np.save('/export/home/tfoulquier/heal_new/MAPCORR.npy',thecorr)

hp.mollview(idata[:,0]*1E6,min=-100,max=100,cmap='jet',hold=False,sub=(2,2,1)   ,nest=True)
hp.mollview(idata[:,1]*1E6,min=-100,max=100,cmap='jet',hold=False,sub=(2,2,2)   ,nest=True)
hp.mollview(thecorr[:,0]*1E6,min=-100,max=100,cmap='jet',hold=False,sub=(2,2,3),nest=True)
hp.mollview(thecorr[:,1]*1E6,min=-100,max=100,cmap='jet',hold=False,sub=(2,2,4),nest=True)
plt.show()

