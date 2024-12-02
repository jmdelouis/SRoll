import healpy as hp
import numpy as np
import foscat.Spline1D as sp1d
import foscat.CircSpline as spC1d
from datetime import datetime, timedelta
import time
import sys
import os
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import cg

def thefunc(x):
    return (-2*x**3+3*x**2)

def compute3deg(x,N):
    if isinstance(x,float):
        ix=int(x*(N-1))
        i=np.zeros([4],dtype='int')
        r=np.zeros([4])
    else:
        ix=(x*(N-1)).astype('int')
        i=np.zeros([4,x.shape[0]],dtype='int')
        r=np.zeros([4,x.shape[0]])
    xx=x*(N-1)-ix
    r[3]=thefunc(xx/2)/2
    r[2]=thefunc(0.5+xx/2)/2
    r[1]=thefunc(1-xx/2)/2
    r[0]=thefunc(0.5-xx/2)/2
    i[3]=ix+3
    i[2]=ix+2
    i[1]=ix+1
    i[0]=ix
    
    if isinstance(x,float):
        if i[0]==0:
            i[0]=1
        if i[1]==0:
            i[1]=1
        if i[2]==N+1:
            i[2]=N
        if i[3]==N+1:
            i[3]=N
        if i[3]==N+2:
            i[3]=N
    else:
        i[0,i[0]==0]=1
        i[1,i[1]==0]=1
        i[2,i[2]==N+1]=N
        i[3,i[3]==N+1]=N
        i[3,i[3]==N+2]=N
        
    i=i-1
    r=r*r
    return i,r/np.sum(r,0)

def compute_3deg_spline(x,N):
    if isinstance(x,float):
        it0_moins=int(x*(N-1))
        iv=np.zeros([2],dtype='int')
        v=np.zeros([2])
        if it0_moins<N-1:
            xx=x*(N-1)-it0_moins
            v[1]=-2*xx**3+3*xx**2
            iv[1]=it0_moins+1
    else:
        it0_moins=(x*(N-1)).astype('int')
        iv=np.zeros([2,x.shape[0]],dtype='int')
        v=np.zeros([2,x.shape[0]])
        idx2=np.where(it0_moins<N-1)[0]
        if idx2.shape[0]>0:
            xx=x[idx2]*(N-1)-it0_moins[idx2]
            v[1,idx2]=-2*xx**3+3*xx**2
            iv[1,idx2]=it0_moins[idx2]+1
        
    xx=it0_moins+1-x*(N-1)
    v[0]=-2*xx**3+3*xx**2
    iv[0]=it0_moins
    return iv,v

def compute_1deg_spline(x,N):
    if isinstance(x,float):
        it0_moins=int(x*(N-1))
        iv=np.zeros([2],dtype='int')
        v=np.zeros([2])
        if it0_moins<N-1:
            v[1]=x*(N-1)-it0_moins
            iv[1]=it0_moins+1
    else:
        it0_moins=(x*(N-1)).astype('int')
        iv=np.zeros([2,x.shape[0]],dtype='int')
        v=np.zeros([2,x.shape[0]])
        idx2=np.where(it0_moins<N-1)[0]
        if idx2.shape[0]>0:
            v[1,idx2]=x[idx2]*(N-1)-it0_moins[idx2]
            iv[1,idx2]=it0_moins[idx2]+1
        
    v[0]=it0_moins+1-x*(N-1)
    iv[0]=it0_moins
    return iv,v

class corrtime:
  def __init__(self,params,rank):
    self.MINTIME  = params['MINTIME']
    self.MAXTIME  = params['MAXTIME']
    self.DELTATIME  = params['DELTATIME']
    self.begin_rg=params['BeginRing']
    self.end_rg=params['EndRing']
    self.nside=params['Nside']
    self.TimeStep=params['TimeStep']
    self.maxspline=np.max([2,int((self.MAXTIME-self.MINTIME)/(24*3600*self.TimeStep))])
    self.maxth=np.sin(0.113)
    
    if rank==0:
      print('Init corrtime')
      sys.stdout.flush()

    self.angref=8/180.*np.pi
    self.solution={}
    self.solution[-1]=0
    self.xref={}
    self.valmin=2*np.pi/180.0
    self.valmax=12*np.pi/180.0
    self.nspline_max=0
    #### New marine 
    self.splineval = {}
    self.doprint=0
    self.rank=rank
    self.c0_min=10.0
    self.c0_max=0.0

    if rank==0:
      print('End init corrtime')
      sys.stdout.flush()

  def eval(self,
           signal,
           hit,
           Inc,
           rg,
           ib,
           hidx,
           externals):
    
    # Incidence 
    ang=np.array(Inc)
    c0 = np.clip((ang-self.valmin)/(self.valmax-self.valmin),0.0,1.0)
            
    # compute time
    t0 = np.array(externals).reshape(len(externals)//2,2)
    # do t0 in [0,1] 
    t0 = np.clip((t0[:,1]-self.MINTIME)/(self.DELTATIME),0.0,1.0)
    #print(t0)
    #Signal 
    Y=np.array(signal)

    # Ring number
    # check the ring distribution
    # Ring number
    rg=np.array(rg)

    t,p=hp.pix2ang(self.nside,hidx[0])
            
    test=0

    if hidx[0] in self.solution:
        if self.solution[hidx[0]] is not None:
            nspline=self.solution[hidx[0]].shape[0]//2
        else:
            nspline=np.max([int(self.maxspline*self.maxth/np.sin(t)),1])
    else:
        nspline=np.max([int(self.maxspline*self.maxth/np.sin(t)),1])
      
    while nspline>1 and test==0:
        
        iv,v=compute3deg(t0,nspline)
        idx=np.arange(t0.shape[0])
        
        if nspline>10:
            ncol=t0.shape[0]
            X=coo_matrix((c0*v[0], (iv[0], idx)), shape=(2*nspline,ncol))+ \
               coo_matrix((v[0], (iv[0]+nspline, idx)), shape=(2*nspline, ncol))+\
               coo_matrix((c0*v[1], (iv[1], idx)), shape=(2*nspline, ncol))+ \
               coo_matrix((v[1], (iv[1]+nspline, idx)), shape=(2*nspline, ncol))+\
               coo_matrix((c0*v[2], (iv[2], idx)), shape=(2*nspline, ncol))+ \
               coo_matrix((v[2], (iv[2]+nspline, idx)), shape=(2*nspline, ncol))+\
               coo_matrix((c0*v[3], (iv[3], idx)), shape=(2*nspline, ncol))+ \
               coo_matrix((v[3], (iv[3]+nspline, idx)), shape=(2*nspline, ncol))
        else:
            X=np.zeros([nspline*2,t0.shape[0]])

            X[iv[0],idx]=v[0]*c0
            X[iv[1],idx]+=v[1]*c0
            X[iv[2],idx]+=v[2]*c0
            X[iv[3],idx]+=v[3]*c0

            X[iv[0]+nspline,idx]=v[0]
            X[iv[1]+nspline,idx]+=v[1]
            X[iv[2]+nspline,idx]+=v[2]
            X[iv[3]+nspline,idx]+=v[3]

        mat=X@X.T
        solution=None
        try:
            if nspline>10:
                solution,e=cg(M,X@Y)
                if e==0:
                    test=1
                else:
                    solution=None
                    nspline-=1
            else:
                cond=np.linalg.cond(mat)
                if cond<1000:
                    solution=np.linalg.solve(mat,X@Y)
                    test=1
                else:  
                    nspline-=1
        except:
            nspline-=1
              
    if nspline==1:
        X=np.ones([2,t0.shape[0]])
        X[0,:]=c0
        mat=X@X.T
        solution=None
        try:
            cond=np.linalg.cond(mat)
            if cond<1000:
                solution=np.linalg.solve(mat,X@Y)
                test=1
        except:
            test=0
          

    if test==0 or solution is None:
        self.solution[hidx[0]]=None
        return [k*0.0 for k in Inc],hit

    corr=solution@X
    del X

    if nspline>1:
        itmp,tmp=compute3deg(0.5,nspline)
        vref=np.zeros([nspline*2])
        for k in range(4):
            vref[itmp[k]]+=0.6*tmp[k]
            vref[nspline+itmp[k]]+=tmp[k]
    else:
        vref=np.zeros([2])
        vref[0]=0.6
        vref[1]=1.0

    xref=solution@vref
    
    self.solution[hidx[0]]=solution
    self.xref[hidx[0]]=xref
      
    corr=corr-xref
    
    if not np.isfinite(np.min(corr)):
        if self.doprint==0:
            self.doprint=1
          
            print('NAN ',nspline,t0,c0,Y)
        self.solution[hidx[0]]=None
        return [k*0.0 for k in Inc],[0.0 for k in hit]

    return [corr[k] for k in range(t0.shape[0])],hit

  def eval_correction(self,
                      inc_ref,
                      rg_ref,
                      hidx):
    if hidx in self.solution:
      if self.solution[hidx] is None:
        return [0.0]
      else:
        nspline=self.solution[hidx].shape[0]//2
        if nspline==1: # no time variation
            return [0.0]
        
        itmp,tmp=compute3deg(rg_ref,nspline)
        vref=np.zeros([nspline*2])
        for k in range(4):
            vref[itmp[k]]+=((inc_ref-self.valmin)/(self.valmax-self.valmin))*tmp[k]
            vref[nspline+itmp[k]]+=tmp[k]
        x=self.solution[hidx]@vref
        return [self.xref[hidx]-x]
    else:
      return [0.0]
    
  def npar_correction(self,
                      hidx):
    if hidx in self.solution:
      if self.solution[hidx] is None:
        return -1
      else:
        return self.solution[hidx].shape[0]
    else:
      return 0
  
class azi:

  def __init__(self,params,rank,begin_rg,end_rg):
    self.MINTIME  = params['MINTIME']
    self.MAXTIME  = params['MAXTIME']
    self.DELTATIME  = params['DELTATIME']
    self.ntime_spline = params['ntime_spline']
    self.histo= params['fhisto']
    self.nspline=params['n_spline']
    self.nbolo=len(params['bolo'])
    self.TimeStep=params['TimeStep']
    self.rank=rank
    self.time_spline_resolution=params['time_spline_resolution']
    if rank==0:
      print('Init azimuth correction')
      sys.stdout.flush()

    ref={}
    idx={}
    myspC1d = spC1d.CircSpline(self.nspline,3)

    for i in range(1000):
      vv2=np.array(myspC1d.calculate(2*np.pi*i/(1000)))
      lidx2=np.where(vv2>0.0)[0]
      idx[i]=[int(lidx2[l]) for l in range(len(lidx2))]
      ref[i]=[vv2[lidx2[l]] for l in range(len(lidx2))]
    self.spline_idx=idx
    self.spline_ref=ref
    n=24*self.TimeStep # ratio between timestep and offset step
    self.kernel=np.ones(n)/n
    if rank==0:
      print('End init azimuth correction')
      sys.stdout.flush()

  def getnumber_of_channels(self):
    return 1

  def eval(self,rg,ib,hidx,idx,inc,externals):
    a1=int(self.time_spline_resolution*self.ntime_spline*(externals[1]-self.MINTIME)/self.DELTATIME)
    if a1>self.time_spline_resolution*self.ntime_spline:
      print('Time problem ',(externals[1]-self.MINTIME)/self.DELTATIME)
      a1=self.time_spline_resolution*self.ntime_spline
    if a1<0:
      a1=0
    a1=self.histo[ib][a1]

    iv,v=compute3deg(a1,self.ntime_spline)
    
    a2=int(1000*(np.fmod(externals[0]+np.pi,2*np.pi)/(2*np.pi)))

    return [int(iv[k]*self.nbolo+ib) for k in range(4)]+[int((k+self.ntime_spline)*self.nbolo+ib) for k in self.spline_idx[a2]], \
      [v[k] for k in range(4)]+self.spline_ref[a2]
      
  def normalize(self,x):
    q=np.array(x[self.ntime_spline*self.nbolo:])
    avv=1E4*np.sum(q.reshape(self.nspline,self.nbolo),0)

    q=np.array(x[:self.ntime_spline*self.nbolo]).reshape(self.ntime_spline,self.nbolo)
    if self.kernel.shape[0]<self.ntime_spline:
        for k in range(self.nbolo):
            q[:,k]=1E2*np.convolve(q[:,k],self.kernel,mode='same')
    else:
        for k in range(self.nbolo):
            q[:,k]=1E4*np.mean(q[:,k])
        
    q=q.reshape(self.nbolo*self.ntime_spline)
    out=[q[i] for i in range(self.nbolo*self.ntime_spline)]+[avv[i%self.nbolo] for i in range(self.nbolo*self.nspline)]
    
    return out


class proj_spline:

  def __init__(self,params,rank,begin_rg,end_rg):
    self.nside=params['Nside']
    self.date_ini = datetime(1970,1,1,0,0,0)
    self.day = params['day']


  def getnumber_of_channels(self):
    return 1

  def eval(self,
           ptg_tuple_2,
           ext_data,
           rg_idx,
           healpix_idx,
           id_bolo,
           hit,
           idx_in_ring,
           signal,
           calib):
    
    return 1,signal, hit, signal, [1.0]

def main():

  
  # DIR DATA
  dir_data = "/home/datawork-cersat-public/cache/project/deepsee2/sroll_input/"

  # Spline 
  SparseFunc = "azi" 
  
  # info_date = np.loadtxt('/home1/datawork/mgallian/SRoll_datarmor/SRoll/sroll4/file_b1.txt',delimiter=',',dtype='str')
  # Nombre de ring à sélectionner, ici 100 ring
  BeginRing = 0 #int(info_date[0])
  EndRing   = 500 #int(info_date[1]) # max ring is 19880
  day = 'A_%dR'%(EndRing) #info_date[2]
  CorrTOD = "corrtime"
  TimeStep= 2
  do_offset=0

  regrid=0

  # nbr spline 
  n_spline=16
  

  # Type de projection (voir class proj)
  projection = "proj_spline"

  #Nom map output
  OMAP='SWIM_2020-2022'

  bolo=['V0_L2S04','V0_L2S06','V0_L2S08','V0_L2S10']

  out_str = '_4dNEPx1_CT'

  # paramètre utilisé pour scat2healpix.py 
  Nside = 512
  RINGSIZE = 100000

  npt_incidence = 128


  # En prendre 1 ring  sur RSTEP : pour run 1, pour debug  plus de 1 
  RSTEP = 1

  # Normaliser le gain 
  NORM_GAIN = 1

  # Plus il est petit, moins on va calculer de pixel. Indice de validité d'un pixel
  seuilcond = 100 # ns = 512

  # Test sur selection 28/48 sur 20/50 : On pourrait mettre 1 comme conditionnement, à tester  

  # Nombre d'itération pour fit le gain 
  NITT = 3
  N_IN_ITT = 500

  # Limite de calcul
  S_IN_ITT = 1E-18

  # Pour chaque détecteur est ce qu'on rajoute une moyenne à nos fit 
  # val_mean = [0.0 for i in range(8)] # 8 splines inci 
  # w_mean = [1E8 for i in range(8)] # 3 zones
  # do_mean = ([1.0] + [0.0]*7)*16 + ([0.0]+[1.0]+[0.0]*6)*16+ ([0.0]*2+[1.0]+[0.0]*5)*16+ ([0.0]*3+[1.0]+[0.0]*4)*16+ ([0.0]*4+[1.0]+[0.0]*3)*16 + \
  #         ([0.0]*5+[1.0]+[0.0]*2)*16 +([0.0]*6+[1.0]+[0.0])*16 +([0.0]*7+[1.0])*16 

  # do_mean  = []
  # val_mean = []
  # w_mean = []

  #  nombre de détecteurs (SWIM : 10 deg, 6 degres ...) ?
  nbolo = len(bolo)

  # bolomask Définir les cartes :  1 ière carte = tous les détecteurs, 2 ième carte 1 ier et 2 ième détecteur 
  bolomask = [1,1,1,1]
  # mapname = ['V0_filt_Full4']+[b for b in bolo]
  mapname = ['2020-2022_Full'] #+[b for b in bolo]

  # Nombre de survey : On peut grouper plusieurs survey ou les prendre 1 par 1
  MAPRINGS = [1]
  # Description des rings
  beg_surv=[BeginRing]
  end_surv=[EndRing]
  name_surv=['']
  

#%% ####################################### OUTPUTS  ############################################

  #dirout='/home/datawork-cersat-public/cache/project/deepsee2/sroll_output/'
  dirout='/home1/scratch/jmdeloui/'

  # Out path
  Out_MAP = [dirout+"/%s_%s%s_J%s"%(OMAP,mapname[i],out_str,day) for i in range(len(bolomask)//nbolo)]
  Out_VEC = [dirout+"/%s_%s%s_J%s"%(OMAP,i,out_str,day) for i in bolo]
  Out_Offset = [dirout+"/%s_%s%s_J%s"%(OMAP,i,out_str,day) for i in bolo]
  Out_Offset_corr = [dirout+"/%s_%s%s_J%s"%(OMAP,i,out_str,day) for i in bolo]
  Out_xi2 = [dirout+"/%s_%s%s_J%s"%(OMAP,i,out_str,day) for i in bolo]
  Out_xi2_corr = [dirout+"/%s_%s%s_J%s"%(OMAP,i,out_str,day) for i in bolo]
  #%% ############################## INPUTS ####################################

  # Signal un par détecteur 
  Signal = ["%s/%s_%s_sig"%(dir_data,OMAP,i) for i in bolo]
  
  # Phase
  External = []
  for i in bolo:
     External =  External + ["%s/%s_%s_ant_azi"%(dir_data,OMAP,i), 
                             "%s/%s_%s_time"%(dir_data,OMAP,i)] 

  # Nombre de fois qu'on est passé dans une cellule healpix pour chaque ring 
  Hit = ["%s/%s_%s_hit"%(dir_data,OMAP,i) for i in bolo]
  
  # Pointage 
  Ptg = ["%s/%s_%s_ptg"%(dir_data,OMAP,i) for i in bolo]
  
  MAXMPIBUFFER = 1024*1024

  #%% ################### Pas utilisé  ######################## 

  # Nan value 
  UNSEEN = -1.6375e30

  # ???
  verbose = 0
  BUILDTF = 0
  SEED= [1,1,1,1]

  # Calibration de départ des détecteurs 
  Calibration = [1.0,1.0,1.0,1.0]

  # Monopole à supprimer des données
  Monop = [0.0,0.0,0.0,0.0]

  # Bruit blanc de chaque détecteurs
  NEP = [3.16,2.449, 2.236, 1.732]

  # if GAINSTEP =0 do not calibrate the data, if GAINTSEP>0 calibrate the data with 1 gain,if GAINSTEP<0 calibrate the dat after -GAINSTEP iterations
  GAINSTEP = 0

  n_bolo =len(bolo)

  ###################  ###################  ###################

  # compute the time scale
  MINTIME=1E30
  MAXTIME=0

  for k in range(len(Hit)):
    r = np.fromfile(External[2*k+1],offset=4*BeginRing*RINGSIZE,count=RINGSIZE,dtype='float32')
    h = np.fromfile(Hit[k],offset=4*BeginRing*RINGSIZE,count=RINGSIZE,dtype='float32')
    a=np.min(r[h>0])-1000.0
    if a<MINTIME:
      MINTIME=a

    r = np.fromfile(External[2*k+1],offset=4*EndRing*RINGSIZE,count=RINGSIZE,dtype='float32')
    h = np.fromfile(Hit[k],offset=4*EndRing*RINGSIZE,count=RINGSIZE,dtype='float32')
    a=np.max(r[h>0])+1000.0
    if a>MAXTIME:
      MAXTIME=a

  DELTATIME=MAXTIME-MINTIME

  ntime_spline=int(DELTATIME/(3600))+1
  time_spline_resolution=16

    
  nhisto=time_spline_resolution*ntime_spline
  histo=np.zeros([nhisto+1])
  
  try:
      from mpi4py import MPI
      comm=MPI.COMM_WORLD
      rank=comm.Get_rank()
      size=comm.Get_size()
      isMPI=1
      
  except:
      rank=0
      size=1
      isMPI=0
      
  print('rank ',rank,size)

  ntimestep=int(DELTATIME/(24*3600))

  fhisto={}
  for k in range(nbolo):
      histo=np.zeros([nhisto+1])
      for rg in range(BeginRing+rank,EndRing+1,size):
          r = np.fromfile(External[2*k+1],offset=4*rg*RINGSIZE,count=RINGSIZE,dtype='float32')
          h = np.fromfile(Hit[k],offset=4*rg*RINGSIZE,count=RINGSIZE,dtype='float32')
          l_time=(nhisto*(r[h>0]-MINTIME)/DELTATIME).astype('int')
          histo=histo+np.bincount(l_time,minlength=nhisto+1)
          
      if isMPI==1:
          histo=comm.allreduce(histo,op=MPI.SUM)

      histo=np.cumsum(histo)
  
      print(rank,k,histo.max(),MINTIME,MAXTIME,DELTATIME,'Number of days ',ntimestep)
    
      fhisto[k]=histo/histo.max()
    
      if rank==0:
          np.save(dirout+"/%s_%s_HISTO.npy"%(OMAP,bolo[k]),histo)
          
  rg_sub_ref=[(k+0.5)/ntimestep for k in range(ntimestep)]
  inc_sub_ref=[8.0/180.0*np.pi for k in range(ntimestep)]
  name_sub=['T%03d'%(k) for k in range(ntimestep)]

  print('reading ok')
  sys.stdout.flush()
  params = vars()
  return params


if __name__ == "__main__":
    main()
