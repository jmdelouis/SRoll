
.. _parameters:

List of parameters
==================


.. list-table:: SRoll parameters
   :widths: 50 50
   
   
   * - NAME
     - DESCRITPION
   * - bolo
     - List of detectors name
   * - nbolo
     - Size of the list 'bolo'
   * - BeginRing
     - First ring to consider (included)
   * - EndRing
     - Last ring to consider (included)
   * - RSTEP
     - Ring step. If set to 1, all the input rings are processed by sroll,if set to 10, 1 ring out of 10 is used, and so on. Useful to quick check some parameter             combinations, as run time is significantly reduced.Use RSTEP=1 for actual productions.
   * - Nside
     - Indice of pixelization of the sphere healpix
   * - D_NOPOL
     -
   * - NORM_GAIN
     - Flag for normalise the gain 
   * - REMOVE_CAL
     - Flag for remove calibration
   * - ADDDIP
     - Flag for add the dipole 
   * - KCMBIN
     -
   * - TESTPOL
     - If set to 0, activates the sroll simulation code (without stim),where colored noise and ADCNL residuals are added to the input signal HPR.If set to 4, only the ADCNL residuals part is added (?). Other values allowed in the code are: 3 (commented), 7 (?). Ignored (forced to -1) if stim_paramfiles is set.
   * - XI2STOP
     - Threshold for the XI2 calcul
   * - seuilcond
     - Threshold used for cond matrix
   * - NITT
     - Number of iterations
   * - FITANGLE
     -
   * - FITPOLEFF
     -
   * - saveCOV
     -
   * - saveCO
     -
   * - delta_psi
     -  Angle in degrees to add to all bolometers psi pointing ({dbpath}/{pixname}_REP6_ptg_TUPLE_2) to simulate an entire focal plane rotation.Omit or set to 0.0 for default/legacy behavior (no focal plane rotation)
   * - verbose
     -
   * - TEMPLATE_NSIDE
     - Indice of pixelisation for the template
   * - BUILDTF
     -
   * - SEED
     - 
   * - GAINSTEP
     - Number of gain values to fit over the whole mission.
       If set to 1, only one gain value will be fitted, thus no variable gains.
       RD12 values are GAINSTEP=128 for 100GHz to 217GHz and GAINSTEP=32 for 353GHz to 857GHz.
   * - NADU
     - 
   * - NADUSTEP
     - Number of setp used in fit_adu ( function no longer use in sroll4)
   * - DOCNN
     - Bolean list of size NITT, run the CNN or not at each iteration
   * - val_mean
     -
   * - w_mean
     - 
   * - do_mean
     - 
   * - CALLCNN
     - Path for CNN
   * - CNN_START
     - Define at wich iteration start to use CNN
   * - CNN_XSIZE
     -
   * - CNN_YSIZE
     -
   * - CNN_RESIDU
     -
   * - CNN_WEIGHTS
     -
   * - N_IN_ITT
     -
   * - INST_CNN
     -
   * - MAP_CNN
     -
   * - Calibration
     - 
   * - NEP
     -
   * - CrossPol
     - Table of Polarization efficiency [%] for each detectors
   * - Monop
     -
   * - OUT_NOPOL
     -
   * - bolomask
     -
   * - beg_surv
     - Define 
   * - end_surv
     -
   * - name_surv
     -
   * - MAPRINGS
     - 
   * - Mask
     -  Binary map that is used to indicate which pixels or regions of a spherical map are considered valid or usable. ( Example 
   * - projectionType
     - Projection use for the output map : example for cosmology 'I,Q,U' for intensity and polarized map and for oceanography 'spline3' to project signal in 0,6,180degrees.
   * - in_template_map
     -
   * - Signal_noPS
     - Input signal 
   * - ADU
     -
   * - Badring
     - List of flagged rings that wont be process
   * - HPR_Calib
     - List of HPRs containing the component that will be subtracted,not fitted, from Signal_noPS.
   * - Hit_noPS
     -
   * - Ptg_noPS
     - Pointing
   * - Theo_HPR
     - 
   * - Theo_MAP
     - 
   * - phase
     - Input phase of the HPR to be projected to maps. Used for advance denoising
   * - rgcnn
     -
   * - bolo_map
     - List of name for the output (len(bolo_map) == len(MAPRING))
   * - Out_MAP
     - Path to save maps ouput
   * - Out_VEC
     - Path to save vectors ouput
