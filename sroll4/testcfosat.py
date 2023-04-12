import healpy as hp

def main():

  nbolo = 4
  bolo=['PRED5_LS04','PRED5_LS06','PRED5_LS08','PRED5_LS10']
  old_bolo=['857-1','857-2','857-3','857-4']

  BeginRing = 0
  EndRing= 499
  RSTEP = 100

  NORM_GAIN = 1
  REMOVE_CAL =1

  delta_psi = 0.0
  D_NOPOL = 0
  KCMBIN = 0
  ADDDIP = 0
  TESTPOL = 1

  CUTRG= 1
  Nside = 64

  XI2STOP = 1.0
  seuilcond = 100
  NITT = 5

  verbose = 0
  TEMPLATE_NSIDE = 128
  BUILDTF = 0
  UNSEEN = hp.UNSEEN




  ################################## #PARAMETRES RELIE AU RESEAU DE NEURONE  ################################## 
  DOCNN = [0,0,0,0,1,1,1,0,0,0,0,0,0]
  CALLCNN = "/home3/homedir7/perso/tfoulqui/workspace/srollexx_work/sroll/py_function_pwst_pol.py"
  CNN_START = NITT
  #PARAMETRES DIMENSION CNN: 12*32*32 carte nside  =12 et 2 canaux Q et U
  CNN_XSIZE = 256
  CNN_YSIZE = 256
  CNN_RESIDU = 0.0
  CNN_WEIGHTS = "/export/home/tfoulquier/WIDXR_files"
  N_IN_ITT = 300

  #Path for netcdf
  INST_CNN = '/home3/homedir7/perso/tfoulqui/workspace/srollexx_work/sroll'
  MAP_CNN =  '/home3/homedir7/perso/tfoulqui/workspace/srollexx_work/sroll'
  ############################################################################################################ 

  Calibration = [1,1,1,1]  
  CrossPol = [0,0,0,0]
  NEP = [1,1,1,1]
  SEED= [1,2,3,4]


  GAINSTEP = 1
  NADU = [16,16,16,16]
  NADUSTEP= [1,1,1,1]

  Monop = [0]
  FSLCOEF = [0.0]
  OUT_NOPOL = [1]
  n_OUT_NOPOL=len(OUT_NOPOL)

  bolomask = [1,1,1,1]


  beg_surv=[0]
  end_surv=[499]
  name_surv=['Full']

  MAPRINGS = [1]
  
  # ############################## INPUTS ####################################


  Mask = "/export/home1/jmdeloui/CFOSAT/CFOSAT_MASK_%d"%(Nside)
  TEMPLATEMAP = "/export/home/tfoulquier/data_sroll/MAP/map_857_2018.float32.bin"


  projectionType = "spline3" # "I"#"I,Q,U"

  in_template_map_I = ["/export/home/tfoulquier/data_sroll/MAP/map_null.float32.bin" for i in range(0,nbolo)]
  in_template_map_Q = ["/export/home/tfoulquier/data_sroll/MAP/map_null.float32.bin"for i in range(0,nbolo)]
  in_template_map_U = ["/export/home/tfoulquier/data_sroll/MAP/map_null.float32.bin"for i in range(0,nbolo)]

  Signal_noPS = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_Signal"%(i) for i in bolo]

  ADU = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_phase"%(i) for i in bolo]

  Badring = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_discarded_rings"%(i) for i in bolo]

  HPR_Calib = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_CalibWW3"%(i) for i in bolo]

  Hit_noPS = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_HitCorr"%(i) for i in bolo]

  Ptg_noPS = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_Ptg"%(i) for i in bolo]
  
  nspline=16
  Theo_noPS = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_%d_S%d"%(i,nspline,j) for j in range(nspline) for i in bolo]
  

  fsl = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_Signal"%(i) for i in bolo]
  phase = ["/export/home1/jmdeloui/DATA4SROLL4/%s_REP6_phregul"%(i) for i in old_bolo]
  rgcnn = ["/export/home/tfoulquier/data_sroll/reduced_MAP/reduced_NULL_rgadutot.int32.bin"for i in range(0,nbolo)]


  val_mean = [0.0 for i in range(4)]

  w_mean = [1E8 for i in range(4)]

  #do_mean = 0 #if = 1 do moyennne des bolometre =0 dans projgrad2 dans troll.c // tableau taille nbolo
  do_mean  = []
  # mean of H0
  # for j in range(4):
  #   do_mean=do_mean+[int((i%4)==j) for i in range(4*nspline)]  

  # ####################################### OUTPUTS  ############################################
  bolo_map = ['CFOSAT']  

  Out_MAP = 
  Out_VEC = 
  Out_Offset = 
  Out_xi2 = 
  Out_xi2_corr = 

  params = vars()

  return params

