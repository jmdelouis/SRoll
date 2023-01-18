exec(open('IMO_4_27.py').read())

BeginRing = SURVBEG['full']
EndRing   = SURVEND['full']

angle = 0.0
tl_template = False
onlyparam=False
doreference=False
dodata=False
iloss=-2
RSTEP = @@RSTEP@@
nocnn=False
cnn2d=False
dopsm=False
dosim=False
bolo=['545-1']
nharm=0
dononoise=False
NITT_CNN = 200
NITT = 10
fitfsl=False
NHIDDEN=4
NCOMP = 32
amp=1.0
REGULAMP = 100.0
VAR2FIT = 1.0
DOFULLYCONNECTED=False
DOFSL=False
TFLEARNCNN2D='FSL_545_INITMODEL'
OUTNAME='857_FSL'

for i in range(len(sys.argv)):
    if sys.argv[i]=='--TFLEARN':
        TFLEARNCNN2D=sys.argv[i+1]
        if rank==0:
            print('TFLEARNCNN2D ',TFLEARNCNN2D)
    if sys.argv[i]=='--amp':
       amp=float(sys.argv[i+1])
    if sys.argv[i]=='--regul':
       REGULAMP=float(sys.argv[i+1])
    if sys.argv[i]=='--var2fit':
       VAR2FIT=float(sys.argv[i+1])
    if sys.argv[i]=='--fitfsl':
       fitfsl=True
    if sys.argv[i]=='--nonoise':
       dononoise=True
    if sys.argv[i]=='--data':
       dodata=True
    if sys.argv[i]=='--dopsm':
       dopsm=True
    if sys.argv[i]=='--harm':
        nharm=int(sys.argv[i+1])
        if rank==0:
            print('nharm',nharm)
    if sys.argv[i]=='--npar':
        NHIDDEN=int(sys.argv[i+1])
        if rank==0:
            print('NHIDDEN',NHIDDEN)
    if sys.argv[i]=='--bolo':
        bolo=[sys.argv[i+1]]
        if rank==0:
            print('bolo',bolo)
    if sys.argv[i]=='--sim':
        dosim=True
        isim=int(sys.argv[i+1])
    if sys.argv[i]=='--isloss':
        iloss=int(sys.argv[i+1])
        if rank==0:
            print('===================================')
            print('==                               ==')
            print('==           ISLOSS = %3d        =='%(iloss))
            print('==                               ==')
            print('===================================')
    if sys.argv[i]=='--dofsl':
        if rank==0:
            print('===================================')
            print('==                               ==')
            print('==           DOFSL               ==')
            print('==                               ==')
            print('===================================')
        DOFSL=True
    if sys.argv[i]=='--onlyparam':
        onlyparam=True
        iloss=-32
    if sys.argv[i]=='--doreference':
        NITT_CNN = 4000
        NITT = 2
        doreference=True
        iloss=-32
    if sys.argv[i]=='--nocnn':
        nocnn=True
    if sys.argv[i]=='--cnn2d':
        cnn2d=True
    if sys.argv[i]=='--rstep':
        RSTEP=int(sys.argv[i+1])
        if rank==0:
            print('===================================')
            print('==                               ==')
            print('==           RSTEP  = %3d        =='%(RSTEP))
            print('==                               ==')
            print('===================================')
    if sys.argv[i]=='--fullyconnected':
        DOFULLYCONNECTED=True
        iloss=-32
    if sys.argv[i]=='--template':
        tl_template=True
        iloss=-32
    if sys.argv[i]=='--angle':
        angle = float(sys.argv[i+1])	

WLoss = 10**(iloss)
NITTCALIB = 5
DORELU=True
rgsize= 27664
Nside = 2048


bolo2       = [BOLOID[b] for b in bolo]
Calibration = [DX11CALIB[b] for b in bolo]
Monop       = [DX11MONOP[b] for b in bolo]
NEP         = [DX11NEP[b] for b in bolo]

if dodata==False and dopsm==False and dosim==False:
    Calibration = [1/amp for b in bolo]
    Monop       = [0.0 for b in bolo]
    NEP         = [1.0 for b in bolo]


if dodata==True or dosim==True:
    ibol=bolo[0]
else:
    ibol='545-1'
   
if cnn2d==True:
    TFLEARN = '@@OUT@@/VECT/%s_offsets_PROD_REP6_RD12_545GHz_CNN2D_DOREF_FSL_%d_-32'%(ibol,NHIDDEN)
else:
    TFLEARN = '@@OUT@@/VECT/%s_offsets_PROD_REP6_RD12_545GHz_CNN_DOREF_FSL_%d_-32'%(ibol,NHIDDEN)
if dosim==True: 
    if cnn2d==True:
        TFLEARN = '@@OUT@@/VECT/%s_offsets_PROD_REP6_RD12_545GHz_SIM000_CNN2D_DOREF_FSL_%d_-32'%(ibol,NHIDDEN)
    else:
        TFLEARN = '@@OUT@@/VECT/%s_offsets_PROD_REP6_RD12_545GHz_SIM000_CNN_DOREF_FSL_%d_-32'%(ibol,NHIDDEN)

NUMBER_TF = @@NUMTF@@

seuilcond      = 1e-06
CrossPol       = [XPOL[b] for b in bolo]
in_polar_fit_Q = ['BINARY:FLOAT32:@@MAP128@@/353_Q.float32.bin' for b in bolo]
in_polar_fit_U = ['BINARY:FLOAT32:@@MAP128@@/353_U.float32.bin' for b in bolo]

Mask = 'BINARY:FLOAT32:@@MASK@@'

Badring = ['BINARY:INT32:@@ROI@@/%s_discarded_rings_dx11'%(b) for b in bolo2]

if (doreference==True) and (not (angle==0)):
    Signal = ['BINARY:FLOAT32:@@PBR@@/fsl_shifted_timeline_%.1e'%(angle) for b in bolo]
else:
    Signal = ['BINARY:FLOAT32:@@PBR@@/%s_REP7_2'%(b) for b in bolo] 

if dosim==True:   
    Signal = ['BINARY:FLOAT32:@@PBRSIM@@/%s_JAN18_stimHPR_%03d.float32.bin'%(b,isim) for b in bolo]
if dodata==True:
    Signal = ['BINARY:FLOAT32:@@PBR@@/%s_REP6_splinefill_JAN20'%(b) for b in bolo]
if dopsm==True:
    Signal = ['BINARY:FLOAT32:@@PBRPSM@@/%s_JAN18sky_nostim_HPR'%(b) for b in bolo]

Hit    = ['BINARY:FLOAT32:@@PBR@@/%s_REP6_hit'%(b) for b in bolo]
Ptg_PH = ['BINARY:FLOAT64:@@PBR@@/%s_REP6_ptg'%(b) for b in bolo]
Ptg_TH = ['BINARY:FLOAT64:@@PBR@@/%s_REP6_ptg_TUPLE_1'%(b) for b in bolo]
Ptg_PSI = ['BINARY:FLOAT64:@@PBR@@/%s_REP6_ptg_TUPLE_2'%(b) for b in bolo]

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

#tmpname=['REP6_R0','REP6_R1','REP6_R2','REP6_R3','REP6_H0','REP6_H1','REP6_H2','REP6_H3']
tmpname=['REP6_H0','REP6_H1','REP6_H2','REP6_H3'] 

if fitfsl==True:			  
    tmpname=['REP6_H0','REP6_H1','REP6_H2','REP6_H3','REP7_2'] 

if tl_template==True:
    Template  = [[[i,'LINEAR',      'BINARY:FLOAT32:@@PBR@@/%s_%s'%(b,i),'ABSOLUT'] for i in tmpname]+
                 [['OOF_NOISE','LINEAR','BINARY:FLOAT32:@@PBR@@/fsl_shifted_timeline_%.1e'%(angle),'ABSOLUT'
                 ]] for b in bolo]
else:
    if nocnn==True:
        if dosim==True and dononoise==True:
            Template  = [[[i,'LINEAR',      'BINARY:FLOAT32:@@PBR@@/%s_%s'%(b,i),'RELATIV'] for i in tmpname] for b in bolo]
        else:
            Template  = [[[i,'LINEAR',      'BINARY:FLOAT32:@@PBR@@/%s_%s'%(b,i),'ABSOLUT'] for i in tmpname] for b in bolo]
    else:
        if cnn2d==True:
            if DOFSL==False:
                TFLEARN='@@PBR@@'+'/'+TFLEARNCNN2D
                par=np.load('%s_%d_param.npy'%(TFLEARN,0))
                nbatch,NHIDDEN=par.shape
                ww=np.load('%s_%d_w_%d.npy'%(TFLEARN,0,0))
                NCOMP=ww.shape[2]//2
                DOFULLYCONNECTED=True
                Template  = [[[i,'LINEAR',      'BINARY:FLOAT32:@@PBR@@/%s_%s'%(b,i),'ABSOLUT'] for i in tmpname]+
                             [['OOF_NOISE','CNN2D',   
                               'BINARY:FLOAT32:@@PBR@@/%s_REP6_phregul'%(b), #x parameter
                               'LINEAR',      # choose how the parameter are used (LINEAR|HIST)
                               [0,2*np.pi],   # data range
                               ['RING',-1,'@@PBR@@/CODEIDX_%s_%s.npy'%(OUTNAME,b)],    # y parameter (RING=one per ring|ONE= one,means only 1D data|HPIX[2|2048]=healpix pixels)
                               NHIDDEN,       # Number of parameters
                               [NCOMP*2,NCOMP*2,NCOMP*2,1],  # number of channels 
                               7,             # kernel size (7 is the only debugged one for circular cnn)
                               'NA',    # data model
                               512,           # Output size from CNN
                               1024,           # Output y size from CNN
	                       -nharm,        # negative nharm low-pass filter (0) no filter
	                       REGULAMP,      # norm of the synthetized correction
                               float(VAR2FIT)
                             ]] for b in bolo]
            else:
                Template  = [[[i,'LINEAR',      'BINARY:FLOAT32:@@PBR@@/%s_%s'%(b,i),'ABSOLUT'] for i in tmpname]+
                             [['OOF_NOISE','CNN2D_FSL',   
                               'BINARY:FLOAT32:@@PBR@@/%s_REP6_phregul'%(b), #x parameter
                               'LINEAR',      # choose how the parameter are used (LINEAR|HIST)
                               [-np.pi,np.pi],   # data range
                               ['RING',0],    # y parameter (RING=one per ring|ONE= one,means only 1D data|HPIX[2|2048]=healpix pixels)
                               NHIDDEN,       # Number of parameters
                               [NHIDDEN,NHIDDEN*2,NHIDDEN*2,1],  # number of channels 
                               7,             # kernel size (7 is the only debugged one for circular cnn)
                               'CIRCULAR',    # data model
                               128,           # Output size from CNN
                               128,           # Output y size from CNN
		               -nharm,        # negative nharm low-pass filter (0) no filter
		               REGULAMP,      # norm of the synthetized correction
		               'BINARY:FLOAT32:@@PBR@@/%s_REP7_2'%(b),
                               '@@PBR@@/rotidx_32_128',
                               'BINARY:FLOAT32:@@PBR@@/imref_%s-32_128.npy'%(b),
                               32,
                               float(VAR2FIT)
                             ]] for b in bolo]
        else:
        
            if DOFSL==True:
                TFLEARN='@@PBR@@'+'/'+TFLEARNCNN2D
                par=np.load('%s_%d_param.npy'%(TFLEARN,0))
                nbatch,NHIDDEN=par.shape
                ww=np.load('%s_%d_w_%d.npy'%(TFLEARN,0,0))
                NCOMP=ww.shape[2]//2
                DOFULLYCONNECTED=True
                Template  = [[[i,'LINEAR',      'BINARY:FLOAT32:@@PBR@@/%s_%s'%(b,i),'RELATIV'] for i in tmpname]+
                             [['OOF_NOISE','CNN1D_FSL',   
                               'BINARY:FLOAT32:@@PBR@@/%s_REP6_phregul'%(b), #x parameter
                               'LINEAR',      # choose how the parameter are used (LINEAR|HIST)
                               [0,2*np.pi],   # data range
                               ['RING',0],    # y parameter (RING=one per ring|ONE= one,means only 1D data|HPIX[2|2048]=healpix pixels)
                               NHIDDEN,       # Number of parameters
                               [NCOMP*2,NCOMP*2,NCOMP*2,1],  # number of channels 
                               4,             # kernel size (7 is the only debugged one for circular cnn)
                               'REGULAR',     # data model
                               128,           # Output size from CNN
		               nharm,         # positive nharm high-pass filter (0) no filter
		               REGULAMP,      # norm of the synthetized correction
                               '@@PBR@@/rot_idx_1D_tl',
                               'BINARY:FLOAT32:@@PBR@@/imref_%s-32_128.npy'%(b),
                               32,
                               float(VAR2FIT)
                             ]] for b in bolo]
            else:
                Template  = [[[i,'LINEAR',      'BINARY:FLOAT32:@@PBR@@/%s_%s'%(b,i),'ABSOLUT'] for i in tmpname]+
                             [['OOF_NOISE','CNN1D',   
                               'BINARY:FLOAT32:@@PBR@@/%s_REP6_phregul'%(b), #x parameter
                               'LINEAR',      # choose how the parameter are used (LINEAR|HIST)
                               [0,2*np.pi],   # data range
                               ['RING',0],    # y parameter (RING=one per ring|ONE= one,means only 1D data|HPIX[2|2048]=healpix pixels)
                               NHIDDEN,       # Number of parameters
                               [NHIDDEN,NHIDDEN*2,NHIDDEN*2,1],  # number of channels 
                               7,             # kernel size (7 is the only debugged one for circular cnn)
                               'CIRCULAR',    # data model
                               128,           # Output size from CNN
		               nharm,         # positive nharm high-pass filter (0) no filter
		               REGULAMP       # norm of the synthetized correction
                             ]] for b in bolo]

if dosim==True and doreference==True and dononoise==False:
    if cnn2d==True:
        if DOFSL==False:
            TFLEARN='/export/home1/jmdeloui/DATA4SROLL4/FSL_857_INITMODEL_4'
            NHIDDEN=4
            NCOMP=16
            DOFULLYCONNECTED=True
            Template  = [[[i,'LINEAR',      'BINARY:FLOAT32:@@PBR@@/%s_%s'%(b,i),'ABSOLUT'] for i in tmpname]+
                         [['OOF_NOISE','CNN2D',   
                           'BINARY:FLOAT32:@@PBR@@/%s_REP6_phregul'%(b), #x parameter
                           'LINEAR',      # choose how the parameter are used (LINEAR|HIST)
                           [0,2*np.pi],   # data range
                           ['RING',-1,'@@PBR@@/CODEIDX_%s_%s.npy'%(OUTNAME,b)],    # y parameter (RING=one per ring|ONE= one,means only 1D data|HPIX[2|2048]=healpix pixels)
                           NHIDDEN,       # Number of parameters
                           [NCOMP*2,NCOMP*2,NCOMP*2,1],  # number of channels 
                           7,             # kernel size (7 is the only debugged one for circular cnn)
                           'NA',    # data model
                           512,           # Output size from CNN
                           1024,           # Output y size from CNN
	                   -nharm,        # negative nharm low-pass filter (0) no filter
	                   REGULAMP,      # norm of the synthetized correction
                           float(VAR2FIT)
                         ]] for b in bolo]
        else:
            Template  = [[[i,'LINEAR',      'BINARY:FLOAT32:@@PBR@@/%s_%s'%(b,i),'ABSOLUT'] for i in tmpname]+
                         [['OOF_NOISE','CNN2D_FSL',   
                           'BINARY:FLOAT32:@@PBR@@/%s_REP6_phregul'%(b), #x parameter
                           'LINEAR',      # choose how the parameter are used (LINEAR|HIST)
                           [-np.pi,np.pi],   # data range
                           ['RING',0],    # y parameter (RING=one per ring|ONE= one,means only 1D data|HPIX[2|2048]=healpix pixels)
                           NHIDDEN,       # Number of parameters
                           [NHIDDEN,NHIDDEN*2,NHIDDEN*2,1],  # number of channels 
                           7,             # kernel size (7 is the only debugged one for circular cnn)
                           'CIRCULAR',    # data model
                           128,           # Output size from CNN
                           128,           # Output y size from CNN
		           -nharm,        # negative nharm low-pass filter (0) no filter
		           REGULAMP,      # norm of the synthetized correction
		           'BINARY:FLOAT32:@@PBR@@/%s_REP7_2'%(b),
                           '@@PBR@@/rotidx_32_128',
                           'BINARY:FLOAT32:@@PBR@@/imref_%s-32_128.npy'%(b),
                           32,
                           float(VAR2FIT)
                         ]] for b in bolo]
    else:
        if DOFSL==True:
            TFLEARN='/export/home1/jmdeloui/DATA4SROLL4/FSL_857_INITMODEL_4'
            NHIDDEN=4
            NCOMP=64
            DOFULLYCONNECTED=True
            Template  = [[[i,'LINEAR',      'BINARY:FLOAT32:@@PBR@@/%s_%s'%(b,i),'ABSOLUT'] for i in tmpname]+
                         [['OOF_NOISE','CNN1D_FSL',   
                           'BINARY:FLOAT32:@@PBR@@/%s_REP6_phregul'%(b), #x parameter
                           'LINEAR',      # choose how the parameter are used (LINEAR|HIST)
                           [0,2*np.pi],   # data range
                           ['RING',0],    # y parameter (RING=one per ring|ONE= one,means only 1D data|HPIX[2|2048]=healpix pixels)
                           NHIDDEN,       # Number of parameters
                           [NCOMP*2,NCOMP*2,NCOMP*2,1],  # number of channels 
                           4,             # kernel size (7 is the only debugged one for circular cnn)
                           'REGULAR',     # data model
                           128,           # Output size from CNN
		           nharm,         # positive nharm high-pass filter (0) no filter
		           REGULAMP,      # norm of the synthetized correction
                           '@@PBR@@/rot_idx_1D_tl',
                           'BINARY:FLOAT32:@@PBR@@/imref_%s-32_128.npy'%(b),
                           32,
                           float(VAR2FIT)     
                         ]] for b in bolo]
        else:
            Template  = [[[i,'LINEAR',      'BINARY:FLOAT32:@@PBR@@/%s_%s'%(b,i),'ABSOLUT'] for i in tmpname]+
                         [['OOF_NOISE','CNN1D',   
                           'BINARY:FLOAT32:@@PBR@@/%s_REP6_phregul'%(b), #x parameter
                           'LINEAR',      # choose how the parameter are used (LINEAR|HIST)
                           [0,2*np.pi],   # data range
                           ['RING',0],    # y parameter (RING=one per ring|ONE= one,means only 1D data|HPIX[2|2048]=healpix pixels)
                           NHIDDEN,       # Number of parameters
                           [NHIDDEN,NHIDDEN*2,NHIDDEN*2,1],  # number of channels 
                           7,             # kernel size (7 is the only debugged one for circular cnn)
                           'CIRCULAR',    # data model
                           128,           # Output size from CNN
		           nharm,         # positive nharm high-pass filter (0) no filter
		           REGULAMP       # norm of the synthetized correction
                         ]] for b in bolo]
            
    TemplateMap = [[['12CO','LINEAR',128,'BINARY:FLOAT32:@@MAP128@@/NEW_cleaned_12CO_K.float32.bin','RELATIV',0.0],
                    ['13CO','LINEAR',128,'BINARY:FLOAT32:@@MAP128@@/NEW_cleaned_13CO_K.float32.bin','RELATIV',0.0]] for b in bolo]
else:
    TemplateMap = [[['12CO','LINEAR',128,'BINARY:FLOAT32:@@MAP128@@/NEW_cleaned_12CO_K.float32.bin','RELATIV',0.0],
                    ['13CO','LINEAR',128,'BINARY:FLOAT32:@@MAP128@@/NEW_cleaned_13CO_K.float32.bin','RELATIV',0.0],
                    ['DUST','LINEAR',128,'BINARY:FLOAT32:@@MAP128@@/Dust_I','RELATIV',0.0],
                    ['POLEFF','POLEFF',128,'BINARY:FLOAT32:@@MAP128@@/353_@POLAR@.float32.bin','ABSOLUT'],
                    ['POLANG','POLANG',128,'BINARY:FLOAT32:@@MAP128@@/353_@POLAR@.float32.bin','ABSOLUT'],
                    ['GAIN','CALIB',2048,'BINARY:FLOAT32:@@MAP2048@@/map_545_2018.float32.bin','RELATIV',0.0]] for b in bolo]
       
if dosim==True and dononoise==True:
    TemplateMap = [[['12CO','LINEAR',128,'BINARY:FLOAT32:@@MAP128@@/NEW_cleaned_12CO_K.float32.bin','RELATIV',0.0],
                    ['13CO','LINEAR',128,'BINARY:FLOAT32:@@MAP128@@/NEW_cleaned_13CO_K.float32.bin','RELATIV',0.0]] for b in bolo]
domodel=False 
#==========================================================================================================================================
# INFORMATION TO REMOVE :
#==========================================================================================================================================
# First column defines the coefficient
# Second colum defines the data
#==========================================================================================================================================
if dodata==True:
    Added  = [[[1.0,'BINARY:FLOAT32:@@PBR@@/%s_REP6_fsl'%(b)]] +
    [[1.0,'BINARY:FLOAT32:@@PBR@@/%s_REP6_diporb_quat'%(b)]] for b in bolo]
else:
    Added  = [[[0.0,'BINARY:FLOAT32:@@PBR@@/%s_REP6_fsl'%(b)]] +
    [[0.0,'BINARY:FLOAT32:@@PBR@@/%s_REP6_diporb_quat'%(b)]] for b in bolo]
if dosim==True and doreference==True:
    Added = [[[-1.0027355308607069/DX11CALIB[b],'BINARY:FLOAT32:@@PBRPSM@@/%s_JAN18sky_nostim_HPR'%(b)]] for b in bolo]
if dosim==True and dononoise==True:
    Added = [[[-1.0027355308607069/DX11CALIB[b],'BINARY:FLOAT32:@@PBRPSM@@/%s_JAN18sky_nostim_HPR'%(b)]] for b in bolo]

OUTCLEANCALIB=False

if cnn2d==True:
    SUF='CNN2D_FSL'
else:
    SUF='CNN_FSL'
    
if onlyparam==True:
    if cnn2d==True:
        SUF='CNN2D_TL_FSL'
    else:
        SUF='CNN_TL_FSL'
if doreference==True:
    if cnn2d==True:
        SUF='CNN2D_DOREF_FSL'
    else:
        SUF='CNN_DOREF_FSL'
if nocnn==True:
    SUF='NORM_FSL'

if tl_template==True:
    SUF='TL_SHIFTED_TEMPLATE_%.1e'%(angle)

if dodata==True:
    SUF='DATA_'+SUF
if dosim==True:
    SUF='SIM%03d_'%(isim)+SUF
if dopsm==True:
    SUF='PSM_'+SUF
if nharm!=0:
    SUF=SUF+'_%d'%(nharm)
if dononoise==True:
    SUF=SUF+'_NN'

if fitfsl==True:
    SUF=SUF+'_FFSL'

MAP = [['@@OUT@@/MAP/%s_RD12_REP6_%s_%d_%d_Full'%(bolo[0],SUF,NHIDDEN,iloss),'NOPOL',[1],SURVBEG['full'],SURVEND['full']],
       ['@@OUT@@/MAP/%s_RD12_REP6_%s_%d_%d_hm1'%(bolo[0],SUF,NHIDDEN,iloss), 'NOPOL',[1],SURVBEG['hm1'],SURVEND['hm1']],
       ['@@OUT@@/MAP/%s_RD12_REP6_%s_%d_%d_hm2'%(bolo[0],SUF,NHIDDEN,iloss), 'NOPOL',[1],SURVBEG['hm2'],SURVEND['hm2']]]

Out_Offset = ['@@OUT@@/VECT/%s_offsets_PROD_REP6_RD12_545GHz_%s_%d_%d'%(b,SUF,NHIDDEN,iloss) for b in bolo]
Out_chi2   = ['@@OUT@@/VECT/%s_chi2_PROD_REP6_RD12_545GHz_%s_%d_%d'%(b,SUF,NHIDDEN,iloss) for b in bolo]

