import numpy as np

class test2D:
  
  def __init__(self,params):
    self.valmin=0.488
    self.valmax=0.838

    nspline=8
    self.SeaMask=np.load('SeaMask.npy')
    
    SPLINE=np.loadtxt('SPLINE%d.txt'%(nspline))
    SPLINEC=np.loadtxt('SPLINE_CIRC_%d.txt'%(nspline))

    ref={}
    idx={}
    
    for i in range(256):
      vv1=SPLINE[i,1:-1]
      lidx1=np.where(vv1>0.0)[0]
      idx[i]={}
      ref[i]={}
      
      for j in range(256):
        vv2=SPLINEC[j,1:-1]
        lidx2=np.where(vv2>0.0)[0]
        idx[i][j]=[int(lidx1[k])+nspline*int(lidx2[l]) for k in range(len(lidx1)) for l in range(len(lidx2))]
        ref[i][j]=[vv1[lidx1[k]]*vv2[lidx2[l]] for k in range(len(lidx1)) for l in range(len(lidx2))]

    self.spline_idx=idx
    self.spline_ref=ref
  
  def eval(self,rg,ib,hidx,inc,externals):

    ioff=int(self.SeaMask[hidx])
    
    a1=int(256*(inc-self.valmin)/(self.valmax-self.valmin))
    a2=int(256*(np.fmod(externals[0]+2*np.pi,2*np.pi)/(2*np.pi)))
    
    return [k+64*ioff for k in self.spline_idx[a1][a2]],self.spline_ref[a1][a2]
  
class test:
  
  def __init__(self,params):
    self.valmin=0.488
    self.valmax=0.838
    self.SeaMask=1.0-np.load('SeaMask.npy') # 1 if sea

    self.nspline=8
    
    SPLINE=np.loadtxt('SPLINE%d.txt'%(self.nspline))

    ref={}
    idx={}
    
    for i in range(256):
      vv1=SPLINE[i,1:-1]
      lidx1=np.where(vv1>0.0)[0]
      
      idx[i]=[int(k) for k in lidx1]
      ref[i]=[vv1[lidx1[k]] for k in range(len(lidx1))]

    self.spline_idx=idx
    self.spline_ref=ref
  
  def eval(self,rg,ib,hidx,inc,externals):
    
    a1=int(256*(inc-self.valmin)/(self.valmax-self.valmin))
    
    return self.spline_idx[a1]+[self.nspline+k for k in range(4)],self.spline_ref[a1]+ \
      [self.SeaMask[hidx]*np.fabs(externals[1]**2*np.cos(2*externals[0])),self.SeaMask[hidx]*np.fabs(externals[1]**2*np.sin(2*externals[0])),
       self.SeaMask[hidx]*np.fabs(externals[2]**2*np.cos(2*externals[0])),self.SeaMask[hidx]*np.fabs(externals[2]**2*np.sin(2*externals[0]))]
  
class proj:
  
  def __init__(self,params):
    self.valmin=0.488
    self.valmax=0.838

  def getnumber_of_channels(self):
    return 1
  
  def eval(self,
           ptg_tuple_2,
           ext_data,
           rg_idx,
           healpix_idx,
           id_bolo):
    
    return [1.0]

def main():
  
  # DIR DATA
  dir_data = "/home1/datawork/mgallian/sroll_data/scat_hpr/"
  
  SparseFunc = "test"
  # Type de projection (voir class proj)
  projection = "proj"
  
  #Nom map output
  OMAP='SCAT3'
  
  #  nombre de détecteurs (SWIM : 10 deg, 6 degres ...) ?
  nbolo = 1
  #bolo=['V1_2']
  bolo=['V4_db_vv']
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
  
  val_mean = [0.0 for i in range(2)] 
  # poids dans la matrice 
  w_mean = [1E8 for i in range(2)]
  # normalize 
  do_mean = [1 for k in range(8)]+[0 for k in range(12)]+[1 for k in range(4)]
  
  
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
  # bolomask = [1,1,1,1,
  #             1,0,0,1,
  #             0,1,1,0]
  
  # Pas utilisé  mais doit être de dim : bolo.size /!\  
  old_bolo=['857-1']


  #%% ############################## INPUTS ####################################

  # Masque pour les continents; peut être 2 masques plus tard pour régler pb continuité comme entre amérique centrale 
  #Mask = "/export/home1/jmdeloui/CFOSAT/CFOSAT_MASK_%d"%(Nside)


  # Signal un par détecteur 
  Signal = ["%s/SCAT_%s_sig"%(dir_data,i) for i in bolo]

  # Phase
  ADU = ["%s/SCAT_%s_phase"%(dir_data,i) for i in bolo]

  # Défini les orbits pas bons 
  # Badring = ["%s/SCAT_%s_discarded_rings"%(dir_data,i) for i in bolo]
  # Badring = ['ALLGOOD']

  # Calibration contient les données wave watch
  #HPR_Calib = ["%s/SCAT_PRED4_LS10_CalibWW3"]
  #HPR_Calib = ['/home1/datawork/jmdeloui/SCAT_V2_ECMWF']

  External = ["%s/SCAT_%s_phase"%(dir_data,i) for i in bolo]+ \
    ["/home1/datawork/jmdeloui/SCAT_%s_ECMWF_U"%(i) for i in bolo]+ \
    ["/home1/datawork/jmdeloui/SCAT_%s_ECMWF_V"%(i) for i in bolo]
  
  # Nombre de fois qu'on est passé dans une cellule healpix pour chaque ring 
  Hit = ["%s/SCAT_%s_hit"%(dir_data,i) for i in bolo]
  
  # Pointage 
  Ptg = ["%s/SCAT_%s_ptg"%(dir_data,i) for i in bolo]
  
  # permet de fit les erreurs en fonction de l'angle alpha
  #HPR = ['/home1/datawork/jmdeloui/SCAT_V2_ECMWF_U','/home1/datawork/jmdeloui/SCAT_V2_ECMWF_V']



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

  # PARAMETRES RELIE AU RESEAU DE NEURONE --> old 
  DOCNN = [0]
  CALLCNN = "/home3/homedir7/perso/tfoulqui/workspace/srollexx_work/sroll/py_function_pwst_pol.py"
  CNN_START = NITT
  #PARAMETRES DIMENSION CNN: 12*32*32 carte nside  =12 et 2 canaux Q et U
  CNN_XSIZE = 256
  CNN_YSIZE = 256
  CNN_RESIDU = 0.0
  CNN_WEIGHTS = "/export/home/tfoulquier/WIDXR_files"

  #Path for netcdf
  INST_CNN = '/home3/homedir7/perso/tfoulqui/workspace/srollexx_work/sroll'
  MAP_CNN =  '/home3/homedir7/perso/tfoulqui/workspace/srollexx_work/sroll'

  D_NOPOL = 1 # Utilisé
  KCMBIN = 0
  ADDDIP = 0
  CUTRG= 1
  DOMAXVRAIE = 1
  XI2STOP = 1.0

  # pas utilisé
  n_bolo =len(bolo)

  # Pas utilisé 
  old_bolo=['857-1']

  # Pas utilisé: pour cosmo 
  Monop = [0]
  FSLCOEF = [0.0]
  OUT_NOPOL = [1]
  n_OUT_NOPOL=len(OUT_NOPOL)

  # Pas utilisé
  TEMPLATEMAP = "/export/home/tfoulquier/data_sroll/MAP/map_857_2018.float32.bin"

  # Pas utilisé
  in_template_map_I = ["/export/home/tfoulquier/data_sroll/MAP/map_null.float32.bin" for i in range(0,nbolo)]
  in_template_map_Q = ["/export/home/tfoulquier/data_sroll/MAP/map_null.float32.bin"for i in range(0,nbolo)]
  in_template_map_U = ["/export/home/tfoulquier/data_sroll/MAP/map_null.float32.bin"for i in range(0,nbolo)]
  # Pas utilisé

  # Pas utilisé --> pour réseau neuronne 
  fsl = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_Signal"%(i) for i in bolo]
  rgcnn = ["/export/home/tfoulquier/data_sroll/reduced_MAP/reduced_NULL_rgadutot.int32.bin"for i in range(0,nbolo)]


  ###################  ###################  ###################

  print('reading ok')
  params = vars()

  return params

