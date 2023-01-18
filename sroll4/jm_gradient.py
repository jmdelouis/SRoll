import numpy as np
import matplotlib.pyplot as plt

nring=100
nside=16
ringsz=12*nside*nside

offset=np.cos(np.arange(nring)/nring*2*np.pi)
offset-=offset.mean()

wnoise=np.tile((1+np.arange(ringsz)/ringsz)**2,nring)*10
ir=np.arange(ringsz*nring,dtype='int')//ringsz
data=offset[ir]+wnoise*np.random.randn(ringsz*nring)
wdata=1/(wnoise**2)

hdata=np.tile(np.arange(ringsz,dtype='int'),nring)


def projdata(data,hdata,wdata,ir):
     imh=np.bincount(hdata,weights=wdata)
     im=np.bincount(hdata,weights=data*wdata)

     res=(data-im[hdata]/imh[hdata])*(1-wdata/imh[hdata])
     b2=np.bincount(ir,weights=res)
     return(b2)

def proj(x,hdata,wdata,ir):
     ldata=x[ir]
     imh=np.bincount(hdata,weights=wdata)
     im=np.bincount(hdata,weights=ldata*wdata)

     res=(ldata-im[hdata]/imh[hdata])*(1-wdata/imh[hdata])
     q2=np.bincount(ir,weights=res)
     q2+=x.sum()
     return(q2)


x0=np.zeros([nring])
b2=projdata(data,hdata,wdata,ir)
ax0=proj(x0,hdata,wdata,ir)
r=b2-ax0
p=r


for itt in range(10):
     delta=(r*r).sum()
     print('Delta %.5g'%(delta))
     ap=proj(p,hdata,wdata,ir)
     alpha=delta/(p*ap).sum()

     x0=x0+alpha*p
     r=r-alpha*ap

     beta=(r*r).sum()/delta

     p=r+beta*p


plt.plot(x0)
plt.plot(offset)
plt.show()


