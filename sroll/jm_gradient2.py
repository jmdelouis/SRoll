import numpy as np
import matplotlib.pyplot as plt
import healpy as hp

nring=1000
nside=64
ringsz=4096

offset=np.cos(np.arange(nring)/nring*1.5*np.pi)
offset-=offset.mean()

wnoise=np.tile((1+np.arange(ringsz)/ringsz)**2,nring)
ir=np.arange(ringsz*nring,dtype='int')//ringsz
template=np.random.randn(ringsz*nring)
ib=ir%4
vv=np.arange(4)-1.5
data=offset[ir]+vv[ib]*template+wnoise*np.random.randn(ringsz*nring)
wdata=1/(wnoise**2)

hdata=np.zeros([ringsz*nring],dtype='int')
for i in range(nring):
     x0=np.cos(i/nring*np.pi*2)*np.sin(i/nring*np.pi*3)
     y0=np.sin(i/nring*np.pi*2)*np.sin(i/nring*np.pi*3)
     z0=np.cos(i/nring*np.pi*3)
     xx=np.random.rand(ringsz,3)*2-1
     xx[:,2]=(-x0*xx[:,0]-y0*xx[:,1])/z0
     xx/=np.sum(xx**2,1).reshape(ringsz,1)

hdata[i*ringsz:(i+1)*ringsz]=hp.vec2pix(nside,xx[:,0],xx[:,1],xx[:,2])


def projdata(data,hdata,wdata,ir):
     imh=np.bincount(hdata,weights=wdata)
     im=np.bincount(hdata,weights=data*wdata)
     res=(data-im[hdata]/imh[hdata])
     imw=np.bincount(hdata,weights=res)
     ima=np.bincount(hdata+ib*12*nside*nside,weights=template)

     res2=(res-wdata/imh[hdata]*imw[hdata])
     res3=(template-wdata/imh[hdata]*ima[hdata+ib*12*nside*nside])
     b2=np.zeros([nring+4])
     b2[0:nring]=np.bincount(ir,weights=res2)
     b2[nring:]=np.bincount(ib,weights=res*res3)
     h2=np.bincount(ir,weights=wdata)
     return(b2,imh,h2)

def proj(x,hdata,wdata,ir,imh):
     ldata=x[ir]+x[nring+ib]*template
     im=np.bincount(hdata,weights=ldata*wdata)
     res=(ldata-im[hdata]/imh[hdata])
     imw=np.bincount(hdata,weights=res)
     ima=np.bincount(hdata+ib*12*nside*nside,weights=template)

     res2=(res-wdata/imh[hdata]*imw[hdata])
     res3=(template-wdata/imh[hdata]*ima[hdata+ib*12*nside*nside])

     q2=np.zeros([nring+4])
     q2[0:nring]=np.bincount(ir,weights=res2)
     q2[nring:]=np.bincount(ib,weights=res*res3)
     q2+=x.sum()
     return(q2)


x0=np.zeros([nring+4])
#x0[0:nring]=offset
#x0[nring:]=vv
b2,imh,h2=projdata(data,hdata,wdata,ir)
ax0=proj(x0,hdata,wdata,ir,imh)

for i in range(nring):
     print(i,b2[i],ax0[i])

r=b2-ax0
p=r

for itt in range(100):
     delta=(r*r).sum()
     print('Delta %.5g'%(delta))
     ap=proj(p,hdata,wdata,ir,imh)
     alpha=delta/(p*ap).sum()

     x0=x0+alpha*p
     r=r-alpha*ap

     beta=(r*r).sum()/delta

     p=r+beta*p

imh=np.bincount(hdata,weights=wdata,minlength=12*nside*nside)
im=np.bincount(hdata,weights=data*wdata,minlength=12*nside*nside)
im2=np.bincount(hdata,weights=(data-x0[ir])*wdata,minlength=12*nside*nside)

plt.figure()
hp.mollview(im/imh,cmap='jet',hold=False,sub=(2,1,1),min=-1,max=1)
hp.mollview(im2/imh,cmap='jet',hold=False,sub=(2,1,2),min=-1,max=1)

plt.figure()
plt.plot(x0[0:nring])
plt.plot(offset)
print(x0[nring:])
print(vv)
plt.show()
