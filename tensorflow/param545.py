BeginRing = 240
EndRing = 26050

NUM_EPOCHS = 25000
EVAL_FREQUENCY = 10
NITT = 7
NITTCALIB = 5 
rgsize=27664
Nside = 2048
RSTEP = 1
verbose = 0

bolo=['545-1','545-2','545-4']
bolo2=['14_545_1','34_545_2','73_545_4']

Calibration = [3.29932298836e-16,3.08792864906e-16,2.64855405172e-16]
Monop       = [3.82437622616e-14,2.57192447573e-14,2.82843952366e-14]
NEP         = [2.06267796946e-16,1.7863247056e-16,1.63152889113e-16]

SEED        = [1234]
NUMBEROFTHREAD = 24

#POLAR TO BE DONE
seuilcond      = 1e-06
CrossPol       = [0.913,0.894,0.89]

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

domodel=False

tmpname_absolut=['REP6_H0','REP6_H1','REP6_H2','REP6_H3','REP7_2','REP6_R0','REP6_R1','REP6_R2','REP6_R3']
tmpname_relativ=[]

Template  = [[[i,'LINEAR','/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/PBR_JMD/%s_%s'%(b,i),'ABSOLUT'] for i in tmpname_absolut]+
             [[i,'LINEAR','/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/PBR_JMD/%s_%s'%(b,i),'RELATIV'] for i in tmpname_relativ]+
             [['ADU','SPLINE1D', '/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/PBR_JMD/%s_REP6_adutot'%(b),32]] for b in bolo]

TemplateMap = [[['12CO','LINEAR',128,'/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/sroll_templates/MAP_JMD_128/NEW_cleaned_12CO_K.float32.bin','RELATIV',0.0],
                ['13CO','LINEAR',128,'/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/sroll_templates/MAP_JMD_128/NEW_cleaned_13CO_K.float32.bin','RELATIV',0.0],
                ['DUST','LINEAR',128,'/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/MAP_JMD_128/Dust_I','RELATIV',0.0],
                ['POLEFF','POLEFF',128,'/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/sroll_templates/MAP_JMD_128/353_@POLAR@.float32.bin','ABSOLUT'],
                ['POLANG','POLANG',128,'/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/sroll_templates/MAP_JMD_128/353_@POLAR@.float32.bin','ABSOLUT'],
                ['GAIN','CALIB',2048,'/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/sroll_templates/MAP_JMD_2048/map_545_2018.float32.bin','RELATIV',0.0]] for b in bolo]
        
#==========================================================================================================================================
# INFORMATION TO REMOVE :
#==========================================================================================================================================
# First column defines the coefficient
# Second colum defines the data
#==========================================================================================================================================
Added  = [[[-1.0,'/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/PBR_JMD/%s_REP6_fsl'%(b)]] +
          [[-1.0, '/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/PBR_JMD/%s_REP6_diporb_quat'%(b)]] for b in bolo]

# Do not remove the calibration template from the map (This should not be the case for the dipole)
OUTCLEANCALIB=False

MAP = [['/scratch/cnt0028/ias1717/SHARED/bware/sroll2r94_jmd/545/MAP/545GHz_RD12_REP6_TENS_ADC2_32','NOPOL',[1,1,1]],
       ['/scratch/cnt0028/ias1717/SHARED/bware/sroll2r94_jmd/545/MAP/545-1_RD12_REP6_TENS_ADC2_32','NOPOL',[1,0,0]],
       ['/scratch/cnt0028/ias1717/SHARED/bware/sroll2r94_jmd/545/MAP/545-2_RD12_REP6_TENS_ADC2_32','NOPOL',[0,1,0]],
       ['/scratch/cnt0028/ias1717/SHARED/bware/sroll2r94_jmd/545/MAP/545-4_RD12_REP6_TENS_ADC2_32','NOPOL',[0,0,1]]]

Out_Offset = ['/scratch/cnt0028/ias1717/SHARED/bware/sroll2r94_jmd/545/VECT/%s_offsets_PROD_REP6_RD12_545GHz_TENS_ADC2'%(b) for b in bolo]
#TO BE DONE
Out_chi2 = ['/scratch/cnt0028/ias1717/SHARED/bware/sroll2r94_jmd/545/VECT/%s_chi2_PROD_REP6_RD12_545GHz_TENS_ADC2'%(b) for b in bolo]

