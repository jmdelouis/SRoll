import numpy as np
import foscat.Spline1D as spl 

class test:
  
  def __init__(self,params,rank,begin_rg,end_rg):
    # read the first ring to compute the minvalue
    self.time_min=params['BeginRing']
    self.time_max=params['EndRing']
    self.nspline=params['nspline']
    self.time_step=128
    
    splinetime=spl.Spline1D(self.nspline)
    
    ref1={}
    idx1={}

    for i in range(self.time_step*self.nspline+1):
      vv1=np.array(splinetime.calculate(i/(self.time_step*self.nspline)))
      lidx1=np.where(vv1>0.0)[0]
      if len(lidx1)>4:
        print(i,len(lidx1))
      idx1[i]=[int(k) for k in lidx1]
      ref1[i]=[vv1[lidx1[k]] for k in range(len(lidx1))]
      
    self.spline_time_idx=idx1
    self.spline_time_ref=ref1
  
  def eval(self,rg,ib,hidx,inc,externals):
      
    a1=int(self.time_step*self.nspline*(rg-self.time_min)/(self.time_max-self.time_min))

    return self.spline_time_idx[a1],self.spline_time_ref[a1]
  
class proj:
  
  def __init__(self,params,rank,begin_rg,end_rg):

    self.rank=rank
    self.l_beg=begin_rg
    self.l_end=end_rg
    self.eta=[(1-k)/(1+k) for k in params['CrossPol']]
    self.Calibration=[1/k for k in params['Calibration']]

  def getnumber_of_channels(self):
    print('get number of chan')
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

    calib=ext_data[id_bolo]

    return signal,hit,calib,[1.0]


def main():
  
  # DIR DATA
  dir_data = "/travail/cwalter/PBR_JMD"

  # Nombre de ring a selectionner, ici 100 ring
  BeginRing = 240
  EndRing = 26050
  
  # Type de projection (voir class proj)
  projection = "proj"
  
  # Class definissant les corrections 
  SparseFunc = "test"
  nspline = 100

  # do the offset
  do_offset=0

  #Nom map output
  OMAP='100-1A'
  
  #  nombre de détecteurs (SWIM : 10 deg, 6 degres ...) ?
  bolo=['100-1a']
  bolo2=['00_100_1a']

  # paramètre utilisé pour scat2healpix.py 
  Nside = 64

  # En prendre 1 ring  sur RSTEP : pour run 1, pour debug  plus de 1 
  RSTEP = 100

  # Sélection de ring par portion
  beg_surv=[BeginRing]
  end_surv=[EndRing]
  name_surv=['Full']

  # Pour chaque détecteur est ce qu'on rajoute une moyenne à nos fit
  
  val_mean = [0.0] 
  # poids dans la matrice 
  w_mean = [1.0]
  # normalize 
  do_mean = [1.0 for k in range(nspline)]
  
  # Nombre de survey : On peut grouper plusieurs survey ou les prendre 1 par 1
  MAPRINGS = [1]
  
  # Normaliser le gain 
  NORM_GAIN = 0

  # Soustraire wave watch 3 à la première itération (NITT) 
  REMOVE_CAL = 1


  # Plus il est petit, moins on va calculer de pixel. Indice de validité d'un pixel
  seuilcond = 100

  # Nombre d'itération pour fit le gain 
  NITT = 5
  
  # Nombre d'itération entre chaque iteration de gain
  N_IN_ITT = 100
  
  # Limite de calcul
  S_IN_ITT = 1E-20
  
  # Nan value 
  UNSEEN = -1.6375e30 

  # 1 seul gain pour toute la mission
  GAINSTEP = 1

  # print some error message
  verbose = 0
  BUILDTF = 0
  SEED= [1]


  # Calibration de départ des détecteurs (coeff + polarisation)
  Calibration = [1.00340513024e-13]  

  # Calibration de départ des détecteurs (polarisation)
  CrossPol = [0.0272]

  # Bruit blanc de chaque détecteurs
  NEP = [2.4767234259e-16]
  
  Monop = [2.2605358109e-15]

  # bolomask défini quel bolometre sont utilisés pour calculer les cartes :  
  bolomask = [1]

  #%% ############################## INPUTS ####################################

  #Mask should be adapted to the local nside
  Mask = '/travail/jdelouis/MASK_SROLL/MASK_%d_100'%(Nside)

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
