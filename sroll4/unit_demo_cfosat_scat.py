import numpy as np
import foscat.Spline1D as spl 
  
class diag_incidence:
  
  def __init__(self,params):

    self.valmin=0.205
    self.valmax=0.945
    self.npt_incidence=params['npt_incidence']
    
  def getnumber_of_index(self):
    return self.npt_incidence
    
  def get_diag_idx(self,rg,ib,hidx,inc,externals):
    a1=int(self.npt_incidence*(inc-self.valmin)/(self.valmax-self.valmin))
    return a1
    
class test:
  
  def __init__(self,params):
    self.valmin=0.205
    self.valmax=0.945
    
    self.nspline=params['nspline']
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
    a1=int(256*self.nspline*(inc-self.valmin)/(self.valmax-self.valmin))
    
    return self.spline_idx[a1],self.spline_ref[a1]
  
class proj:
  
  def __init__(self,params):
    self.valmin=0.205
    self.valmax=0.945
    self.npt_incidence=params['npt_incidence']

    self.std=1/np.fromfile('NOISEMOD',dtype=float)

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
           signal):
    
    a1=int(self.npt_incidence*(ptg_tuple_2-self.valmin)/(self.valmax-self.valmin))
    if a1<0 or a1>=self.npt_incidence:
      hit=0
      
    hit=hit*self.std[a1]
    
    return signal,hit,[1.0]

def main():
  
  # DIR DATA
  dir_data = "/home1/datawork/mgallian/sroll_data/scat_hpr_v3.3.1"

  # function describing a diag
  DiagFunc = "diag_incidence"
  npt_incidence=128
  
  # function describing the parameters to fit for instrumental correction
  SparseFunc = "test"
  nspline=8
  
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

  # Nombre de ring à sélectionner, ici 100 ring
  BeginRing = 0
  EndRing= 3

  # En prendre 1 ring  sur RSTEP : pour run 1, pour debug  plus de 1 
  RSTEP = 1

  # Sélection de ring par portion
  beg_surv=[BeginRing]
  end_surv=[EndRing]
  name_surv=['Full']

  # Pour chaque détecteur est ce qu'on rajoute une moyenne à nos fit
  do_mean = []
  
  val_mean = [0.0 for i in range(1)] 
  # poids dans la matrice 
  w_mean = [1E8 for i in range(1)]
  # normalize 
  do_mean = [1 for k in range(nspline)]
  
  
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
  N_IN_ITT = 1000
  
  # Nan value 
  UNSEEN = -1.6375e30 
  #hp.UNSEEN

  # 1 seul gain pour toute la mission
  GAINSTEP = 0
  # On sait pas ce que c'est mais 1 si 1 capteur 
  NADU = [1]
  NADUSTEP= [1]

  # ???
  verbose = 0
  BUILDTF = 0
  SEED= [1]


  # Calibration de départ des détecteurs (coeff + polarisation)
  Calibration = [1.0]  

  # Calibration de départ des détecteurs (polarisation)
  CrossPol = [0.0]

  # Bruit blanc de chaque détecteurs
  NEP = [1]
  
  # bolomask Définir les cartes :  1 ière carte = tous les détecteurs, 2 ième carte 1 ier et 2 ième détecteur 
  bolomask = [1]

  #%% ############################## INPUTS ####################################

  # Masque pour les continents; peut être 2 masques plus tard pour régler pb continuité comme entre amérique centrale 
  #Mask = "/export/home1/jmdeloui/CFOSAT/CFOSAT_MASK_%d"%(Nside)


  # Signal un par détecteur 
  Signal = ["%s/SCAT_%s_sig"%(dir_data,i) for i in bolo]

  #External = ["%s/SCAT_%s_phase"%(dir_data,i) for i in bolo]
  
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
  Out_Offset = [dirout+"/%s%s"%(OMAP,i) for i in bolo]
  Out_Offset_corr = [dirout+"/%s%s"%(OMAP,i) for i in bolo]
  Out_xi2 = [dirout+"/%s%s"%(OMAP,i) for i in bolo]
  Out_xi2_corr = [dirout+"/%s%s"%(OMAP,i) for i in bolo]


  #%% ################### Pas utilisé  ########################
  
  D_NOPOL = 1 # Utilisé
  KCMBIN = 0
  ADDDIP = 0
  CUTRG= 1
  DOMAXVRAIE = 1
  XI2STOP = 1.0

  # Pas utilisé: pour cosmo 
  Monop = [0]
  FSLCOEF = [0.0]
  OUT_NOPOL = [1]
  n_OUT_NOPOL=len(OUT_NOPOL)

  ###################  ###################  ###################

  print('reading ok')
  
  params = vars()

  a=proj(params)
  
  return params

if __name__ == "__main__":
    main()
