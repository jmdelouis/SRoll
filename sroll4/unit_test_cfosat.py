import healpy as hp

def main():

 
  nbolo = 4
  bolo=['PRED5_LS04','PRED5_LS06','PRED5_LS08','PRED5_LS10']
  
  old_bolo=['857-1','857-2','857-3','857-4']
  BeginRing = 0
  EndRing= 499
  RSTEP = 100

  beg_surv=[240,240,13144,240,5720,11194,16691,21720]
  end_surv=[26050,13145,26050,5721,11195,16692,21721,27005]
  name_surv=['Full','hm1','hm2','s1','s2','s3','s4','s5']
  MAPRINGS = [3,1,1,1,1]

  n_bolo =len(bolo)

  NORM_GAIN = 1
  REMOVE_CAL = 1

  delta_psi = 0.0
  D_NOPOL = 0
  KCMBIN = 0
  ADDDIP = 0
  TESTPOL = 1

  CUTRG= 1
  DOMAXVRAIE = 1
  Nside = 64
  
  XI2STOP = 1.0
  seuilcond = 100.0
  NITT = 10
  verbose = 0
  TEMPLATE_NSIDE = 128
  BUILDTF = 0
  UNSEEN = hp.UNSEEN
  
  SEED= [1,2,3,4]
  GAINSTEP = 1
  NADU = [16,16,16,16]
  NADUSTEP= [1,1,1,1]

  #PARAMETRES RELIE AU RESEAU DE NEURONE
  DOCNN = [0,0,0,0,1,1,1,0,0,0,0,0,0]
  #DOCNN = [0,0,0,0,0,0,0,0,0,0,0,0,0]  

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

  #Calibration = [3.30076826046e-16,3.55811287601e-16,3.18681631353e-16,2.219187708e-16]
  
  Calibration = [1.0,1,1,1]  

  CrossPol = [0.0,0.0,0.0,0.0]

 
  NEP = [1.0,1.0, 1.0,1.0]
  Monop = [0,0,0,0]
  FSLCOEF = [0.0,0.0,0.0,0.0]
  OUT_NOPOL = [1,1,1,1,1]
  n_OUT_NOPOL=len(OUT_NOPOL)
   # bolomask
  bolomask = [1,1,1,1,
              1,0,0,1,
              0,1,1,0]
  
  MAPRINGS = [19,1,1]

  ##**##
  ############################################## INPUTS ############################################################
  
  Mask = "/export/home1/jmdeloui/CFOSAT/CFOSAT_MASK_%d"%(Nside)
  TEMPLATEMAP = "/export/home/tfoulquier/data_sroll/MAP/map_857_2018.float32.bin"
  
  projectionType = "spline3" # "I"#"I,Q,U"

  in_template_map_I = ["/export/home/tfoulquier/data_sroll/MAP/map_null.float32.bin" for i in range(0,nbolo)]
  
  in_template_map_Q = ["/export/home/tfoulquier/data_sroll/MAP/map_null.float32.bin"for i in range(0,nbolo)]
  
  in_template_map_U = ["/export/home/tfoulquier/data_sroll/MAP/map_null.float32.bin"for i in range(0,nbolo)]
  
  #Signal_noPS = ["/export/home/tfoulquier/data_sroll/MAP/%s_REP6_THEO_SIMU"%(i) for i in bolo]
  #Signal_noPS = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_Signal"%(i) for i in bolo]
  Signal_noPS = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_SignalCorr"%(i) for i in bolo]
  OMAP='3'
  
  #Signal_noPS = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_CalibWW3"%(i) for i in bolo]
  #OMAP='WW3_3'
  

  
  ADU = ["/export/home/tfoulquier/data_sroll/MAP/%s_REP6_adutot"%(i) for i in old_bolo]
  
  Badring = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_discarded_rings"%(i) for i in bolo]

  #DipOrb_noPS = ["/export/home1/jmdeloui/DATA4SROLL4/%s_REP6_diporb_quat"%(i) for i in old_bolo]
  HPR_Calib = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_CalibWW3"%(i) for i in bolo]


    
  Hit_noPS = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_HitCorr"%(i) for i in bolo]
  
  Ptg_noPS = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_Ptg"%(i) for i in bolo]
  
  nspline=16
  Theo_noPS = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_%d_S%d"%(i,nspline,j) for j in range(nspline) for i in bolo]
  
  #do_mean = 0 #if = 1 do moyennne des bolometre =0 dans projgrad2 dans troll.c // tableau taille nbolo
  val_mean = [0.0 for i in range(4)]
  
  w_mean = [1E8 for i in range(4)]
  
  do_mean  = []
  # mean of H0
  for j in range(4):
    do_mean=do_mean+[int((i%4)==j) for i in range(4*nspline)]   
  
  fsl = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_Signal"%(i) for i in bolo]
  phase = ["/export/home1/jmdeloui/DATA4SROLL4/%s_REP6_phregul"%(i) for i in old_bolo]
  rgcnn = ["/export/home/tfoulquier/data_sroll/reduced_MAP/reduced_NULL_rgadutot.int32.bin"for i in range(0,nbolo)]



  ####################################### OUTPUTS  ##################################################################
  bolo_map = ['CFOSAT','LS04_10','LS06_08']
  MAP = ["/export/home/tfoulquier/MAP/CFOSAT_%s%s"%(OMAP,i) for i in bolo_map]
  Out_VEC = ["/export/home/tfoulquier/VEC/CFOSAT_%s%s_offset"%(OMAP,i) for i in bolo]
  Out_Offset = ["/export/home/tfoulquier/VEC/CFOSAT_%s%s_offset"%(OMAP,i) for i in bolo]
  Out_Offset_corr = ["/export/home/tfoulquier/VEC/CFOSAT_%s%s_offset_corr"%(OMAP,i) for i in bolo]
  Out_xi2 = ["/export/home/tfoulquier/VEC/CFOSAT_%s%s_xi2"%(OMAP,i) for i in bolo]
  Out_xi2_corr = ["/export/home/tfoulquier/VEC/CFOSAT_%s%s_xi2_corr"%(OMAP,i) for i in bolo]


  params = vars()

  return params




