import numpy as np
import foscat.Spline1D as spl 
  
class diag_incidence:
  
  def __init__(self,params,rank,beginring,endring):

    self.valmin=0.205
    self.valmax=0.945
    self.npt_incidence=params['npt_incidence']
    
  def getnumber_of_index(self):
    return self.npt_incidence
    
  def get_diag_idx(self,rg,ib,hidx,inc,externals):
    a1=int(self.npt_incidence*(inc-self.valmin)/(self.valmax-self.valmin))
    return a1


class test:
  
  def __init__(self,params,rank,beginring,endring):
    # read the first ring to compute the minvalue
    self.time_step=256
    tmp=np.fromfile(params['External'][0],offset=params['BeginRing']*4*3000000,count=3000000,dtype='float32')
    htmp=np.fromfile(params['Hit'][0],offset=params['BeginRing']*4*3000000,count=3000000,dtype='float32')

    if rank==0:
      print('Begin Time ',tmp[htmp>0].min())
    self.timemin=tmp[htmp>0].min()
    del(tmp)
    del(htmp)
    tmp=np.fromfile(params['External'][0],offset=params['EndRing']*4*3000000,count=3000000,dtype='float32')
    htmp=np.fromfile(params['Hit'][0],offset=params['EndRing']*4*3000000,count=3000000,dtype='float32')
    
    if rank==0:
      print('End Time ',tmp[htmp>0].max())
    self.timemax=tmp[htmp>0].max()
    del(tmp)
    del(htmp)
    
    self.valmin=0.205
    self.valmax=0.945
    self.nspline=params['nspline']
    
    self.nsplinetime=params['nspline_time']
    if rank==0:
      print('Use %d knots spline to compute the time variation'%(self.nsplinetime))
    splinetime=spl.Spline1D(self.nsplinetime)
    
    ref1={}
    idx1={}

    for i in range(self.time_step*self.nsplinetime+1):
      vv1=np.array(splinetime.calculate(i/(self.time_step*self.nsplinetime)))
      lidx1=np.where(vv1>0.0)[0]
      if len(lidx1)>4:
        print(i,len(lidx1))
      # AJOUTE NSPLINE POUR EVITER DE CONFONDRE LES INDEXS
      idx1[i]=[int(k+self.nspline) for k in lidx1]
      ref1[i]=[vv1[lidx1[k]] for k in range(len(lidx1))]
      
    self.spline_time_idx=idx1
    self.spline_time_ref=ref1
    
    myspl=spl.Spline1D(self.nspline)

    ref={}
    idx={}

    for i in range(256*self.nspline):
      vv1=np.array(myspl.calculate(i/(256*self.nspline)))
      lidx1=np.where(vv1>0.0)[0]
      idx[i]=[int(k) for k in lidx1]
      ref[i]=[vv1[lidx1[k]] for k in range(len(lidx1))]

    self.spline_idx=idx
    self.spline_ref=ref
  
  def eval(self,rg,ib,hidx,inc,externals):
    
    if inc>self.valmax:
      print(inc)
    if inc<self.valmin:
      print(inc)
      
    if externals[0]>self.timemax:
      print('PROBLEM IN TIME MAX ',externals[0])
    if externals[0]<self.timemin:
      print('PROBLEM IN TIME MIN ',externals[0])
      
    a1=int(256*self.nspline*(inc-self.valmin)/(self.valmax-self.valmin))
    
    atime=int(self.time_step*self.nsplinetime*(externals[0]-self.timemin)/(self.timemax-self.timemin))

    if atime<0 or atime>self.nsplinetime*self.time_step:
      print('PROBLEM IN TIME SCALE')
      print(atime,self.nsplinetime*self.time_step,externals[0],self.timemax)
      exit(0)
    
    return self.spline_idx[a1]+self.spline_time_idx[atime],self.spline_ref[a1]+self.spline_time_ref[atime]
  
class proj:
  
  def __init__(self,params,rank,beginring,endring):
    self.valmin=0.205
    self.valmax=0.945
    self.npt_incidence=params['npt_incidence']

    self.std=1/np.fromfile('NOISEMOD',dtype=float)

  def getnumber_of_channels(self):
    return 2
  
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
      
    a1=int(self.npt_incidence*(ptg_tuple_2-self.valmin)/(self.valmax-self.valmin))
    if a1<0 or a1>=self.npt_incidence:
      hit=0
      
    hit=hit/self.std[a1]
    
    return signal,hit,calib,[1.0,np.cos(ptg_tuple_2)**2]

def main():
  
  # DIR DATA
  dir_data = "/home1/datawork/mgallian/sroll_data/scat_hpr_v3.3.1"

  # Nombre de ring à sélectionner, ici 100 ring
  BeginRing = 0
  EndRing = 30
  
  # function describing a diag
  DiagFunc = "diag_incidence"
  npt_incidence=128
  
  # function describing the parameters to fit for instrumental correction
  SparseFunc = "test"
  nspline=8
  do_offset=0
  nspline_time=(1+EndRing- BeginRing)*2
  
  # Type de projection (voir class proj)
  projection = "proj"
  
  #Nom map output
  OMAP='SCAT3'
  
  #  nombre de détecteurs (SWIM : 10 deg, 6 degres ...) ?
  nbolo = 1
  #bolo=['V1_2']
  bolo=['V8_db_vv']
  inci_str=''

  # paramètre utilisé pour scat2healpix.py 
  Nside = 512

  # Resolution healpix Du template WW3 
  TEMPLATE_NSIDE = Nside

  # En prendre 1 ring  sur RSTEP : pour run 1, pour debug  plus de 1 
  RSTEP = 1

  # Sélection de ring par portion
  beg_surv=[BeginRing]
  end_surv=[EndRing]
  name_surv=['Full']

  # Pour chaque détecteur est ce qu'on rajoute une moyenne à nos fit
  
  val_mean = [0.0 for i in range(2)] 
  # poids dans la matrice 
  w_mean = [1.0,1.0]
  # normalize 
  do_mean = [1.0 for k in range(nspline)]+[0.0 for k in range(nspline_time)]+[0.0 for k in range(nspline)]+[1.0 for k in range(nspline_time)]
  
  # Nombre de survey : On peut grouper plusieurs survey ou les prendre 1 par 1
  MAPRINGS = [1]
  
  # Normaliser le gain 
  NORM_GAIN = 1

  # Soustraire wave watch 3 à la première itération (NITT) 
  REMOVE_CAL = 1

  # Rotation angle
  delta_psi = 0.0

  # Plus il est petit, moins on va calculer de pixel. Indice de validité d'un pixel
  seuilcond = 10.0

  # Nombre d'itération pour fit le gain 
  NITT = 1
  
  # Nombre d'itération entre chaque iteration de gain
  N_IN_ITT = 1000
  
  # Limite de calcul
  S_IN_ITT = 1E-30
  
  # Nan value 
  UNSEEN = -1.6375e30 
  #hp.UNSEEN

  # 1 seul gain pour toute la mission
  GAINSTEP = 0

  # ???
  verbose = 0
  BUILDTF = 0
  SEED= [1]


  # Calibration de départ des détecteurs (coeff + polarisation)
  Calibration = [1.0]  

  # Bruit blanc de chaque détecteurs
  NEP = [1.0]
  
  # bolomask Définir les cartes :  1 ière carte = tous les détecteurs, 2 ième carte 1 ier et 2 ième détecteur 
  bolomask = [1]

  #%% ############################## INPUTS ####################################

  # Masque pour les continents; peut être 2 masques plus tard pour régler pb continuité comme entre amérique centrale 
  #Mask = "/export/home1/jmdeloui/CFOSAT/CFOSAT_MASK_%d"%(Nside)


  # Signal un par détecteur 
  Signal = ["%s/SCAT_%s_sig"%(dir_data,i) for i in bolo]

  External = ["%s/SCAT_%s_time"%(dir_data,i) for i in bolo]
  
  # Nombre de fois qu'on est passé dans une cellule healpix pour chaque ring 
  Hit = ["%s/SCAT_%s_hit"%(dir_data,i) for i in bolo]
  
  # Pointage 
  Ptg = ["%s/SCAT_%s_ptg"%(dir_data,i) for i in bolo]



  #%% ####################################### OUTPUTS  ############################################
  bolo_map = ['SCAT']  

  dirout='/home1/datawork/jmdeloui/scat_maps/'

  # Out path
  Out_MAP = [dirout+"/%s_%s%s"%(OMAP,i,inci_str) for i in bolo]
  Out_VEC = [dirout+"/%s%s"%(OMAP,i) for i in bolo]
  
  # Pas utilisé: pour cosmo 
  Monop = [0]

  ###################  ###################  ###################

  
  params = vars()

  return params

if __name__ == "__main__":
    main()
