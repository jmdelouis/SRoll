import numpy as np

class proj:
  
  def __init__(self,params,rank,begin_rg,end_rg):

    self.rank=rank
    self.l_beg=begin_rg
    self.l_end=end_rg
    self.eta=[(1-k)/(1+k) for k in params['CrossPol']]

  def getnumber_of_channels(self):
    print('get number of chan')
    return 3
  
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
    
    calib=ext_data[id_bolo]

    return signal,hit,calib,[1.0,
                             self.eta[id_bolo]*np.cos(2*ptg_tuple_2),
                             self.eta[id_bolo]*np.sin(2*ptg_tuple_2)]
             
def main():
  
  # DIR DATA
  dir_data = "/travail/cwalter/PBR_JMD"

  # Nombre de ring a selectionner, ici 100 ring
  BeginRing = 240
  EndRing = 26050
  
  # Type de projection (voir class proj)
  projection = "proj"
  
  # do the offset
  do_offset=1

  #Nom map output
  OMAP='100ds1'
  
  #  nombre de détecteurs (SWIM : 10 deg, 6 degres ...) ?
  bolo=['100-1a','100-1b','100-4a','100-4b']
  bolo2=['00_100_1a','01_100_1b','80_100_4a','81_100_4b']

  # Resolution carte finale
  Nside = 32

  # En prendre 1 ring  sur RSTEP : pour run 1, pour debug  plus de 1 
  RSTEP = 400

  # Sélection de ring par portion pour carte en survey a la fin
  beg_surv=[BeginRing]
  end_surv=[EndRing]
  name_surv=['Full']

  # Pour chaque parametres sparse ajoute-t-on une contrainet sur la moyenne
  val_mean = [] 
  # poids dans la matrice 
  w_mean = []
  # normalize 
  do_mean = []
  
  # Nombre de survey dans cahcune des cartes produites. 
  # La definition des surveys est pris dans les tableaux beg_surv,end_surv  etc.
  MAPRINGS = [1]
  
  # Normaliser le gain entre detecteur (cas de calibration absolue inconnue) 
  NORM_GAIN = 0

  # Soustrait le calibrateur dans la carte produite à la fin
  REMOVE_CAL = 1


  # Plus il est petit, moins on va calculer de pixel. Indice de validité d'un pixel
  seuilcond = 100

  # Nombre d'itération pour fit le gain 
  NITT = 4
  
  # Nombre d'itération entre chaque iteration de gain
  N_IN_ITT = 500
  
  # Limite de calcul
  S_IN_ITT = 1E-24
  
  # Nan value 
  UNSEEN = -1.6375e30 

  # 1 seul gain pour toute la mission
  GAINSTEP = 1

  # print some error message
  verbose = 0
  BUILDTF = 0
  SEED= [1]


  # Calibration de départ des détecteurs (coeff + polarisation)
  Calibration = [1.00340513024e-13,1.23240311873e-13,1.46826316998e-13,1.17469574177e-13]  

  # Calibration de départ des détecteurs (polarisation)
  CrossPol = [0.0272,0.0293,0.0219,0.0402]

  # Bruit blanc de chaque détecteurs
  NEP = [2.4767234259e-16,2.56729781043e-16,2.26421345033e-16,2.45551518841e-16]
  
  Monop = [2.2605358109e-15,6.33346357328e-15,6.29208169223e-16,3.36970326506e-15]

  # bolomask Définir les cartes :  1 ière carte = tous les détecteurs, 2 ième carte 1 ier et 2 ième détecteur 
  bolomask = [1,1,1,1]

  #%% ############################## INPUTS ####################################

  #Mask should be adapted to the local nside
  Mask = '/travail/jdelouis/MASK_SROLL/MAKS_%d_100'%(Nside)

  #badring definition
  Badring= ['/travail/cwalter/calROIs/%s_discarded_rings_dx11'%(k) for k in bolo2]

  # Signal un par détecteur 
  Signal = ["%s/%s_REP6"%(dir_data,i) for i in bolo]
  
  # Nombre de fois qu'on est passé dans une cellule healpix pour chaque ring 
  Hit = ["%s/%s_REP6_hit"%(dir_data,i) for i in bolo]
  
  # Pointage 
  Ptg = ["%s/%s_REP6_ptg"%(dir_data,i) for i in bolo]

  # External 
  External =["%s/%s_dipHFI17_quad_hprbin"%(dir_data,i) for i in bolo]


  #%% ####################################### OUTPUTS  ############################################

  dirout='/travail/jdelouis/sroll4_map/'

  # Out path
  Out_MAP = [dirout+"/%s"%(OMAP)]
  Out_VEC = [dirout+"/%s_VEC_%s"%(OMAP,i) for i in bolo]


  ###################  ###################  ###################

  
  params = vars()

  return params

if __name__ == "__main__":
    main()
