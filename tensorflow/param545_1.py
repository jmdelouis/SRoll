BeginRing = 240
EndRing = 26050

NUM_EPOCHS = 10000
EVAL_FREQUENCY = 10
NITT = 1

XIMAGE_SIZE = 128
NDCONV = 3
NHIDDEN  = 8
DEEPNESS = 8
KERNELSZ = 5
NUM_DCONV_CHAN=[ DEEPNESS, DEEPNESS, DEEPNESS, 1]

rgsize=27664
Nside = 2048
RSTEP = 1
verbose = 0

bolo=['545-1']
bolo2=['14_545_1']

Calibration = [3.29932298836e-16]
Monop       = [3.82437622616e-14]
NEP         = [2.06267796946e-16]

SEED        = [1234]

NUMBEROFTHREAD = 4


seuilcond      = 1e-06
CrossPol       = [0.913]
in_polar_fit_Q = ['/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/sroll_templates/MAP_JMD_128/353_Q.float32.bin' for b in bolo]
in_polar_fit_U = ['/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/sroll_templates/MAP_JMD_128/353_U.float32.bin' for b in bolo]

Mask = '/scratch/cnt0028/ias1717/SHARED//RD12_data/dmc/MISS03/DATA/MASK_2048_GALACTIC/mask_RD_545'

Badring = ['/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/calROIs/%s_discarded_rings_dx11'%(b) for b in bolo2]

Signal = ['/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/PBR_JMD/%s_REP6'%(b) for b in bolo]
Hit    = ['/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/PBR_JMD/%s_REP6_hit'%(b) for b in bolo]
Ptg_PH = ['/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/PBR_JMD/%s_REP6_ptg'%(b) for b in bolo]
Ptg_TH = ['/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/PBR_JMD/%s_REP6_ptg_TUPLE_1'%(b) for b in bolo]
Ptg_PSI = ['/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/PBR_JMD/%s_REP6_ptg_TUPLE_2'%(b) for b in bolo]

#==========================================================================================================================================
# SYSTEMATIC TO FIT :
#==========================================================================================================================================
# Can be described by a Timeline information in Template or by a map information in TemplateMap
#==========================================================================================================================================
# Template/TemplateMap:
# * First column = 'LINEAR|CALIB|ADU' : LINEAR fit and remove the template, CALIB fit and use as calibration pattern, ADU fit using the ADCNL model
# * Second column = The data 
# * Third column =  'NA' if nothing to do, List of 3 value (e.g. [0.0,0,len(bolo)]) indicates that fitted amplitude mean from detector 0 to len(bolo)
#   is equal to 0.0 
#==========================================================================================================================================

tmpname=['REP6_R0','REP6_R1','REP6_R2','REP6_R3','REP6_H0','REP6_H1','REP6_H2','REP6_H3','REP7_2']
#tmpname=['REP6_H0','REP6_H1','REP6_H2','REP6_H3'] #,'REP7_2']
#tmpname=['REP6_R0','REP6_R1','REP6_R2','REP6_R3','REP6_H0','REP6_H1','REP6_H2','REP6_H3']
#tmpname=['REP6_H0']

Template  = [[[i,'LINEAR','/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/PBR_JMD/%s_%s'%(b,i),'ABSOLUT'] for i in tmpname]+
             [['ADU','SPLINE1D',   '/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/PBR_JMD/%s_REP6_adutot'%(b),32]] for b in bolo]

TemplateMap = [[['12CO','LINEAR',128,'/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/sroll_templates/MAP_JMD_128/NEW_cleaned_12CO_K.float32.bin','RELATIV',0.0],
                ['13CO','LINEAR',128,'/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/sroll_templates/MAP_JMD_128/NEW_cleaned_13CO_K.float32.bin','RELATIV',0.0],
                ['DUST','LINEAR',128,'/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/MAP_JMD_128/Dust_I','RELATIV',0.0],
                ['POLEFF','POLEFF',128,'/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/sroll_templates/MAP_JMD_128/353_@POLAR@.float32.bin','ABSOLUT'],
                ['POLANG','POLANG',128,'/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/sroll_templates/MAP_JMD_128/353_@POLAR@.float32.bin','ABSOLUT'],
                ['GAIN','CALIB',2048,'/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/sroll_templates/MAP_JMD_2048/map_545_2018.float32.bin','RELATIV',0.0]] for b in bolo]
       
domodel=False 
#==========================================================================================================================================
# INFORMATION TO REMOVE :
#==========================================================================================================================================
# First column defines the coefficient
# Second colum defines the data
#==========================================================================================================================================
Added  = [[[-0.0,'/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/PBR_JMD/%s_REP6_fsl'%(b)]] +
          [[-1.0, '/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/PBR_JMD/%s_REP6_diporb_quat'%(b)]] for b in bolo]

OUTCLEANCALIB=False

MAP = [['/scratch/cnt0028/ias1717/SHARED/bware/sroll2r94_jmd//545/MAP/545-1_RD12_REP6_TENS','NOPOL',[1]]]

Out_Offset = ['/scratch/cnt0028/ias1717/SHARED/bware/sroll2r94_jmd//545/VECT/%s_offsets_PROD_REP6_RD12_545GHz_TENS'%(b) for b in bolo]
Out_chi2 = ['/scratch/cnt0028/ias1717/SHARED/bware/sroll2r94_jmd//545/VECT/%s_chi2_PROD_REP6_RD12_545GHz_TENS'%(b) for b in bolo]

