#Script that create NET_INFO files with given parameters

# IMPORT
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import sys
#import tensorflow as tf
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()
import time
import os


###INIT
nbolo=4
xsize=256
ysize=256
nout = 32


############################ CNN_FSL ############################
XIMAGE_SIZE=256
YIMAGE_SIZE=256
NCOMP= 4
SCALE=2
NDCONV=2 
DEEPNESS=5
KERNELSZ=5
NHIDDEN = NCOMP
DOFULLYCONNECTED=False
XNHIDDEN=XIMAGE_SIZE // (SCALE**NDCONV)
YNHIDDEN=YIMAGE_SIZE // (SCALE**NDCONV)
NTOTHIDDEN=NHIDDEN*(XNHIDDEN*YNHIDDEN)
NUM_DCONV_CHAN={}

for i in range(NDCONV):
    NUM_DCONV_CHAN[i]=int(DEEPNESS)
NUM_DCONV_CHAN[0] =NHIDDEN
NUM_DCONV_CHAN[NDCONV]=1
INPARAM=4


############################ CNN_MAP ############################
mappar=tf.constant(np.load('/export/home/jmdeloui/heal_cnn/MAPCNN_PAR.npy'))
print("mappar ="+str(mappar))
mapchantab=np.load('/export/home/jmdeloui/heal_cnn/MAPCNN_ARCHI.npy')
print(mapchantab)
mapdepth=mapchantab.shape[0]-1
mapchan=mapchantab[0]
wwmap={}
wbmap={}
for i in range(mapdepth):
    wwmap[i]=tf.Variable(np.load('/export/home/jmdeloui/heal_cnn/MAPCNN_W%d.npy'%(i)))
    wbmap[i]=tf.Variable(np.load('/export/home/jmdeloui/heal_cnn/MAPCNN_B%d.npy'%(i)))


##PRINT for Debug
print("mapchantab "+str(mapchantab)+"\n")
print("NUM_DCONV_CHAN"+str(NUM_DCONV_CHAN)+"\n")

print("mapdepth "+str(mapdepth))
print("wwmap "+str(wwmap))
print("wbmap "+str(wbmap))



######################## SAVE  ########################
CNN_names=['FSL','MAP']

for name in CNN_names:
  if name == "FSL" :
    #Init INFO_CNN
    INFO=np.zeros([11],dtype='float')
    INFO[0]=nbolo    
    INFO[1]=XIMAGE_SIZE
    INFO[2]=YIMAGE_SIZE
    INFO[3]=NDCONV
    INFO[4]=int(DOFULLYCONNECTED==True)
    INFO[5]=DEEPNESS
    INFO[6]=NCOMP
    INFO[7]=SCALE
    INFO[8]=KERNELSZ
    INFO[9]=SCALE
    INFO[10] = KERNELSZ

  if name == "MAP":
    #Init INFO_CNN
    INFO=np.zeros([11],dtype='float')
    INFO[0]= 2    
    INFO[1]= 12 * nout * nout
    INFO[2]=1 
    INFO[3]=mapdepth
    INFO[4]=int(DOFULLYCONNECTED==True)
    INFO[5]=mapchan
    INFO[6]=mapchan
    INFO[7]=4 
    INFO[8]=4
    INFO[9] = 1  
    INFO[10] = 1  

  #Create net_info
  np.save('CNN_%s_NET_INFO.npy'%(name),INFO)



