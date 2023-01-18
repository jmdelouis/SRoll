import numpy as np
import matplotlib.pyplot as plt


corr=np.fromfile('/scratch/cnt0028/ias1717/SHARED/bware/sroll2r94_jmd/VEC/JAN18_ALL_DATA_R10_100-1a_offset_0_2_CORRNL',dtype='float').reshape(2,31,320,320)

for i in range(15):
    plt.subplot(5,6,1+i)
    plt.contourf(((corr[0,i,:,:]/corr[1,i,:,:])*1E6).transpose(),vmin=-10,vmax=10)

corr=np.fromfile('/scratch/cnt0028/ias1717/SHARED/bware/sroll2r94_jmd/VEC/JAN18_ALL_DATA_R1_100-1a_offset_0_2_CORRNL',dtype='float').reshape(2,31,320,320)

for i in range(15):
    plt.subplot(5,6,16+i)
    plt.contourf(((corr[0,i,:,:]/corr[1,i,:,:])*1E6).transpose(),vmin=-10,vmax=10)
plt.show()

