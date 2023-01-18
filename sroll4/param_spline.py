import healpy as hp

def main():

 
  nbolo = 4
  #bolo=['PRED2_LS04','PRED2_LS06','PRED2_LS08','PRED2_LS10']
  bolo=['PRED_LS04','PRED_LS06','PRED_LS08','PRED_LS10']
  old_bolo=['857-1','857-2','857-3','857-4']
  BeginRing = 0
  EndRing= 542
  #EndRing= 99
  RSTEP = 1

  n_bolo =len(bolo)

  AVFFF = 857
  AVGFREEFREE = 0
  AVFSYNC = 857
  AVGSYNCHRO  = 0
  AVGDUST100  = 0
  AVFDUST = 857
  AVGDUST = 0
  AVFCO = 857
  AVG12CO = 0
  AVG13CO = 0

  CALCODUST= 0

  CUTRG= 1
  DODIPCAL = 1
  DOGAINDIP = 1
  DOMAXVRAIE = 1
  D_NOPOL = 0
  Nside = 64
  REMDIP = 1
  REMHDIP = 0
  ADDDIP = 0
  KCMBIN = 0
  SAVEINTMAP = 0
  TESTPOL = 1
  XI2STOP = 1.0
  seuilcond = 1E30
  NITT = 3
  FITANGLE = 0	
  FITPOLEFF = 0	
  saveCOV = 0
  saveCO = 0
  delta_psi = 0.0
  verbose = 0
  TEMPLATE_NSIDE = 128
  #POSSIBLE BUG
  #BUILDTF     = 4
  BUILDTF = 0
  UNSEEN = hp.UNSEEN
  
  number_of_SEED = 1
  SEED= 1234

  number_of_ADDCO13= 0

  GAINSTEP = 1

  number_of_NADU = 4
  # 857GHz
  NADU = [16,16,16,16]

  number_of_NADUSTEP = 4
  # 857GHz
  NADUSTEP= [1,1,1,1]


  do_meanBolo0 = 0 #if = 1 do moyennne des bolomettre =0 dans projgrad2 dans troll.c // tableau taille nbolo


  #do_mean_detector - nb 
  #do_mean_listofpix 
  

  
  #PARAMETRES RELIE AU RESEAU DE NEURONE

  DOCNN = [0,0,0,0,1,1,1,0,0,0,0,0,0]
  #DOCNN = [0,0,0,0,0,0,0,0,0,0,0,0,0]  
  number_of_DOCNN = len(DOCNN)

  CALLCNN = "/home3/homedir7/perso/tfoulqui/workspace/srollexx_work/sroll/py_function_pwst_pol.py"

  CNN_START = NITT +1
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
  
  Calibration = [1,1,1,1]  
  number_of_Calibration =len(Calibration)    

  CrossPol = [0.0,0.0,0.0,0.0]
  number_of_CrossPol =len(CrossPol)

 
  NEP = [1.0,1.0, 1.0,1.0]
  number_of_NEP =len(NEP)

  Monop = [0,0,0,0]
  number_of_Monop =len(Monop)
  # 857GHz

  FSLCOEF = [0.0,0.0,0.0,0.0]
  number_of_FSLCOEF=len(FSLCOEF)

  OUT_NOPOL = [1,1,1,1,1]
  n_OUT_NOPOL=len(OUT_NOPOL)
  number_of_OUT_NOPOL=len(OUT_NOPOL)

  # bolomask for 857ghz map
  bolomask = [1,1,1,1,
              1,0,0,0,
              0,1,0,0,
              0,0,1,0,
              0,0,0,1]
  number_of_bolomask=len(bolomask)

  MAPRINGS = [19,1,1,1,1]
  number_of_MAPRINGS = len(MAPRINGS)

  """
  #limite en ring - begin / end ring 
  beginMask = [-1,-1,1200,-1,-1] # Full//HM1/HM2/ODD/EVEN
  endMask = [-1,1200,-1,-1,-1]# Full//HM1/HM2/ODD/EVEN
  #ODDEVEN - 0 pas de flitrage - 1 ring pair / 2 - ring impairs
  ODDEVEN_Mask = [0,0,0,1,2]
  """
  
  ### INPUTS
  Mask = "/export/home1/jmdeloui/CFOSAT/CFOSAT_MASK"
  TEMPLATEMAP = "/export/home/tfoulquier/data_sroll/MAP/map_857_2018.float32.bin"
  #TEMPLATEMAP = "/export/home/jmdeloui/reduced_FAKE_Q"

  #Theo_Dust_Q       = /export/home/tfoulquier/data_sroll/MAP/Dust_New_Q 
  #Theo_Dust_U       = /export/home/tfoulquier/data_sroll/MAP/Dust_New_U 
  #Theo_TDust_I =
  #Theo_TDust_Q =
  #Theo_TDust_U =

  #Theo_FREEFREE =
  #in_synchro_map_I =
  #in_synchro_map_Q =
  #in_synchro_map_U =

  projectionType = "spline3" #"Q,U" #"cfosatmap" #"I"#"I,Q,U"

  in_template_map_I = ["/export/home/tfoulquier/data_sroll/MAP/map_857_2018.float32.bin" for i in range(0,nbolo)]
  number_of_in_template_map_I=len(in_template_map_I)

  in_template_map_Q = ["/export/home/tfoulquier/data_sroll/MAP/map_null.float32.bin"for i in range(0,nbolo)]
  number_of_in_template_map_Q = len(in_template_map_Q)

  in_template_map_U = ["/export/home/tfoulquier/data_sroll/MAP/map_null.float32.bin"for i in range(0,nbolo)]
  number_of_in_template_map_U=len(in_template_map_U)

  in_polar_fit_Q = ["/export/home/tfoulquier/data_sroll/MAP/353_Q.float32.bin"for i in range(0,nbolo)]
  number_of_in_polar_fit_Q = len(in_polar_fit_Q)
  
  in_polar_fit_U = ["/export/home/tfoulquier/data_sroll/MAP/353_U.float32.bin"for i in range(0,nbolo)]
  number_of_in_polar_fit_U = len(in_polar_fit_U)

  #Signal_noPS = ["/export/home/tfoulquier/data_sroll/MAP/%s_REP6_THEO_SIMU"%(i) for i in bolo]
  #Signal_noPS = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_Signal"%(i) for i in bolo]
  Signal_noPS = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_SignalCorr"%(i) for i in bolo]
  #Signal_noPS = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_CalibWW3"%(i) for i in bolo]
  number_of_Signal_noPS = len(Signal_noPS)

  ADU = ["/export/home/tfoulquier/data_sroll/MAP/%s_REP6_adutot"%(i) for i in old_bolo]
  number_of_ADU = len(ADU)

  Badring = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_discarded_rings"%(i) for i in bolo]
  number_of_Badring = len(Badring)

  #DipOrb_noPS = ["/export/home1/jmdeloui/DATA4SROLL4/%s_REP6_diporb_quat"%(i) for i in old_bolo]
  DipOrb_noPS = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_CalibWW3"%(i) for i in bolo]
  number_of_DipOrb_noPS =len(DipOrb_noPS) #-- Calib-temp


  Hit_noPS = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_HitCorr"%(i) for i in bolo]
  number_of_Hit_noPS=len(Hit_noPS)

  Ptg_noPS = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_Ptg"%(i) for i in bolo]
  number_of_Ptg_noPS = len(Ptg_noPS)
  
  #Theo_noPS = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_H%d"%(i,j) for j in range(8) for i in bolo]
  #Theo_noPS = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_8_S%d"%(i,j) for j in range(8) for i in bolo]
  Theo_noPS = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_16_S%d"%(i,j) for j in range(16) for i in bolo]
  number_of_Theo_noPS = len(Theo_noPS)
  
  fsl = ["/export/home1/jmdeloui/CFOSAT/CFOSAT_%s_Signal"%(i) for i in bolo]
  number_of_fsl= len(fsl)


  phase = ["/export/home1/jmdeloui/DATA4SROLL4/%s_REP6_phregul"%(i) for i in old_bolo]
  number_of_phase = len(phase)


  rgcnn = ["/export/home/tfoulquier/data_sroll/reduced_MAP/reduced_NULL_rgadutot.int32.bin"for i in range(0,nbolo)]
  number_of_rgcnn = len(rgcnn)

  ### OUTPUTS
  bolo_map = ['CFOSAT','LS04','LS06','LS08','LS10']
  number_of_bolo_map = len(bolo_map)

  MAP = ["/export/home/tfoulquier/MAP/CFOSAT_%s"%(i) for i in bolo_map]
  number_of_MAP = len(MAP)

  Out_Offset = ["/export/home/tfoulquier/VEC/CFOSAT_%s_offset"%(i) for i in bolo]
  number_of_Out_Offset = len(Out_Offset)


  Out_Offset_corr = ["/export/home/tfoulquier/VEC/CFOSAT_%s_offset_corr"%(i) for i in bolo]
  number_of_Out_Offset_corr = len(Out_Offset_corr)

  Out_xi2 = ["/export/home/tfoulquier/VEC/CFOSAT_%s_xi2"%(i) for i in bolo]
  number_of_Out_xi2 = len(Out_xi2)

  Out_xi2_corr = ["/export/home/tfoulquier/VEC/CFOSAT_%s_xi2_corr"%(i) for i in bolo]
  number_of_Out_xi2_corr = len(Out_xi2_corr)


  params  =vars()

  return params




