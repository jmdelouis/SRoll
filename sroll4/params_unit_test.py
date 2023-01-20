import healpy as hp

def main():

 
  nbolo = 4
  bolo=['857-1','857-2','857-3','857-4']
  bolo2= ['25_857_1','35_857_2','65_857_3','74_857_4']
  BeginRing = 2
  EndRing= 260
  RSTEP = 100

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
  Nside = 2048
  REMDIP = 1
  REMHDIP = 1
  ADDDIP = 0
  KCMBIN = 0
  SAVEINTMAP = 0
  TESTPOL = 1
  XI2STOP = 1.0
  seuilcond = 1e-06
  NITT = 10
  FITANGLE = 0	
  FITPOLEFF = 0	
  saveCOV = 0	
  saveCO = 0
  delta_psi = 0.0
  verbose = 0
  TEMPLATE_NSIDE = 1024
  #POSSIBLE BUG
  #BUILDTF     = 4
  BUILDTF = 0
  #UNSEEN = hp.UNSEEN
  
  number_of_SEED = 1
  SEED= [1,2,3,4]

  number_of_ADDCO13= 0

  GAINSTEP = 1

  number_of_NADU = 4
  # 857GHz
  NADU = [16,16,16,16]

  number_of_NADUSTEP = 4
  # 857GHz
  NADUSTEP= [1,1,1,1]

  #PARAMETRES RELIE AU RESEAU DE NEURONE

  DOCNN = [0,0,0,0,0,0,0]
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

  CrossPol = [0.887,0.883,0.856,0.897]
  number_of_CrossPol =len(CrossPol)

 
  NEP = [2.17663807827,2.36372259024, 2.08546999122,2.01804359351]
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
  bolomask = [1,1,1,1,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]
  number_of_bolomask=len(bolomask)

  MAPRINGS = [19,1,1,1,1]
  number_of_MAPRINGS = len(MAPRINGS)

  """
  #limite en ring - begin / end ring 
  beginMask = [-1,-1,1200,-1,-1] # Full//HM1/HM2/ODD/EVEN
  endMask = [-1,1200,-1,-1,-1]# Full//HM1/HM2/ODD/EVEN
  #ODDEVEN - 0 pas de flitrage - 1 ring pari / 2 - ring impairs
  ODDEVEN_Mask = [0,0,0,1,2]
  """
  
  ### INPUTS
  Mask = "/export/home/tfoulquier/data_sroll/MAP/mask_RD_857"
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

  projectionType = "I" #"I"#"I,Q,U"

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
  Signal_noPS = ["/export/home/jmdeloui/reduced_%s_REP6_FAKE_QU"%(i) for i in bolo]
  number_of_Signal_noPS = len(Signal_noPS)

  ADU = ["/export/home/tfoulquier/data_sroll/reduced_MAP/reduced_%s_REP6_adutot"%(i) for i in bolo]
  number_of_ADU = len(ADU)

  Badring = ["/export/home/tfoulquier/data_sroll/reduced_MAP/reduced_%s_discarded_rings_dx11"%(i) for i in bolo2]
  number_of_Badring = len(Badring)

  DipOrb_noPS = ["/export/home/tfoulquier/data_sroll/reduced_MAP/reduced_%s_REP6_diporb_quat"%(i) for i in bolo]
  number_of_DipOrb_noPS =len(DipOrb_noPS) #-- Calib-temp


  Hit_noPS = ["/export/home/tfoulquier/data_sroll/reduced_MAP/reduced_%s_REP6_hit"%(i) for i in bolo]
  number_of_Hit_noPS=len(Hit_noPS)

  Ptg_noPS = ["/export/home/tfoulquier/data_sroll/reduced_MAP/reduced_%s_REP6_ptg"%(i) for i in bolo]
  number_of_Ptg_noPS = len(Ptg_noPS)

  Theo_noPS = ["/export/home/tfoulquier/data_sroll/reduced_MAP/reduced_%s_REP6_H0"%(i) for i in bolo]
  number_of_Theo_noPS = len(Theo_noPS)


  fsl = ["/export/home/tfoulquier/data_sroll/reduced_MAP/reduced_%s_REP6_fsl"%(i) for i in bolo]
  number_of_fsl= len(fsl)


  phase = ["/export/home/tfoulquier/data_sroll/reduced_MAP/reduced_%s_REP6_phregul"%(i) for i in bolo]
  number_of_phase = len(phase)


  rgcnn = ["/export/home/tfoulquier/data_sroll/reduced_MAP/reduced_NULL_rgadutot.int32.bin"for i in range(0,nbolo)]
  number_of_rgcnn = len(rgcnn)

  ### OUTPUTS
  bolo_map = ['857ghz','857-1','857-2','857-3','857-4']
  number_of_bolo_map = len(bolo_map)

  MAP = ["/export/home/tfoulquier/MAP/sroll2_debug_857ghz_%s"%(i) for i in bolo_map]
  number_of_MAP = len(MAP)

  Out_Offset = ["/export/home/tfoulquier/VEC/sroll2_debug_857ghz_%s_offset"%(i) for i in bolo]
  number_of_Out_Offset = len(Out_Offset)


  Out_Offset_corr = ["/export/home/tfoulquier/VEC/sroll2_debug_857ghz_%s_offset_corr"%(i) for i in bolo]
  number_of_Out_Offset_corr = len(Out_Offset_corr)

  Out_xi2 = ["/export/home/tfoulquier/VEC/sroll2_debug_857ghz_%s_xi2"%(i) for i in bolo]
  number_of_Out_xi2 = len(Out_xi2)

  Out_xi2_corr = ["/export/home/tfoulquier/VEC/sroll2_debug_857ghz_%s_xi2_corr"%(i) for i in bolo]
  number_of_Out_xi2_corr = len(Out_xi2_corr)


  params  =vars()

  return params




