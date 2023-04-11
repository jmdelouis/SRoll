import healpy as hp

def main():

 
  nbolo = 4
  bolo=['857-1','857-2','857-3','857-4']
  bolo2= ['25_857_1','35_857_2','65_857_3','74_857_4'] # ID for badring
  BeginRing = 0
  EndRing= 499
  RSTEP = 100

  beg_surv=[240,240,13144,240,5720,11194,16691,21720]
  end_surv=[26050,13145,26050,5721,11195,16692,21721,27005]
  name_surv=['Full','hm1','hm2','s1','s2','s3','s4','s5']
  
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

  ##################### CNN ################################
  #PARAMETRES RELIE AU RESEAU DE NEURONE
  DOCNN = [0,0,0,1,1,1]
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

  #############################################################
  
  
  #Calibration = [3.30076826046e-16,3.55811287601e-16,3.18681631353e-16,2.219187708e-16]
  
  Calibration = [1.0,1.0,1.0,1.0]  
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
  
  MAPRINGS = [1]

  ##**##
  ############################################## INPUTS ############################################################
  
  Mask = "/home1/scratch/jmdeloui/DATA4SROLL4/mask_RD_857"
  TEMPLATEMAP = "/export/home/tfoulquier/data_sroll/MAP/map_857_2018.float32.bin"
  
  projectionType = "I" #"I,Q,U"

  in_template_map_I = ["/export/home/tfoulquier/data_sroll/MAP/map_857_2018.float32.bin" for i in range(nbolo)]
  
  in_template_map_Q = ["/export/home/tfoulquier/data_sroll/MAP/map_null.float32.bin"for i in range(0,nbolo)]
  
  in_template_map_U = ["/export/home/tfoulquier/data_sroll/MAP/map_null.float32.bin"for i in range(0,nbolo)]
  
  Signal_noPS =  ["/export/home1/jmdeloui/DATA4SROLL4/%s_REP7_2"%(i) for i in bolo]
    
  ADU = ["/export/home1/jmdeloui/DATA4SROLL4/%s_REP6_phregul"%(i) for i in bolo]
  
  Badring = ["/export/home1/jmdeloui/DATA4SROLL4/%s_discarded_rings_dx11"%(i) for i in bolo2]

  HPR_Calib = ["/export/home1/jmdeloui/DATA4SROLL4/%s_REP6_diporb_quat"%(i) for i in bolo]  

  Hit_noPS = ["/export/home1/jmdeloui/DATA4SROLL4/%s_REP6_hit"%(i) for i in bolo]
  
  Ptg_noPS = ["/export/home1/jmdeloui/DATA4SROLL4/%s_REP6_ptg"%(i) for i in bolo]

  OMAP = '3'
  #do_mean = 0 #if = 1 do moyennne des bolometre =0 dans projgrad2 dans troll.c // tableau taille nbolo
  val_mean = [0.0 for i in range(4)]
  
  w_mean = [1E8 for i in range(4)]
  
  nspline = 16
  do_mean  = []
  # mean of H0
  for j in range(4):
    do_mean=do_mean+[int((i%4)==j) for i in range(4*nspline)]   


  Theo_MAP = ['/export/home1/jmdeloui/DATA4SROLL4/NEW_cleaned_12CO.float32.bin',
              '/export/home1/jmdeloui/DATA4SROLL4/FREE_FREE_mod.float32.bin',
              '/export/home1/jmdeloui/DATA4SROLL4/Dust_New_I']

  # remove for polar
  Theo_HPR = ["/export/home1/jmdeloui/DATA4SROLL4/%s_REP6_PDUST"%(i) for i in bolo]
  Theo_HPR = Theo_HPR + ["/export/home1/jmdeloui/DATA4SROLL4/%s_REP6_ODUST"%(i) for i in bolo]
  Theo_HPR = Theo_HPR +["/export/home1/jmdeloui/DATA4SROLL4/%s_REP6_GRAD"%(i) for i in bolo]

  phase = ["/export/home1/jmdeloui/DATA4SROLL4/%s_REP6_phregul"%(i) for i in bolo]
  rgcnn = ["/export/home1/jmdeloui/DATA4SROLL4/%s_REP6_rgadutot.int32.bin"%(i) for i in bolo]




  ####################################### OUTPUTS  ##################################################################
  bolo_map = ['857GHz']

  Out_MAP = 
  Out_VEC = 
  Out_Offset = 
  Out_Offset_corr = 
  Out_xi2 = 
  Out_xi2_corr = 


  params = vars()

  return params




