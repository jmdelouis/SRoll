
.. _parameters:

List of parameters
==================


.. list-table:: SRoll parameters
   :widths: 5 5 5
   :header-rows: 1
   
   * - NAME
     - DESCRITPION
     - TYPE
   * - bolo
     - List of detectors name
     - <class 'list'>
   * - bolo2
     -
     - <class 'list'>
   * - nbolo
     - Size of the list 'bolo'
     - <class 'int'>
   * - BeginRing
     - First ring to consider (included)
     - <class 'int'>
   * - EndRing
     - Last ring to consider (included)
     - <class 'int'>
   * - RSTEP
     - Ring step. If set to 1, all the input rings are processed by sroll,if set to 10, 1 ring out of 10 is used, and so on. Useful to quick check some parameter             combinations, as run time is significantly reduced.Use RSTEP=1 for actual productions.
     - <class 'int'>
   * - D_NOPOL
     -
     - <class 'int'>
   * - NORM_GAIN
     -
     - <class 'int'>
   * - REMOVE_CAL
     -
     - <class 'int'>
   * - ADDDIP
     -
     - <class 'int'>
   * - KCMBIN
     -
     - <class 'int'>
   * - SAVEINTMAP
     -
     - <class 'int'>
   * - TESTPOL
     - If set to 0, activates the sroll simulation code (without stim),where colored noise and ADCNL residuals are added to the input signal HPR.If set to 4, only the ADCNL residuals part is added (?). Other values allowed in the code are: 3 (commented), 7 (?). Ignored (forced to -1) if stim_paramfiles is set.
     - <class 'int'>
   * - XI2STOP
     -
     - <class 'float'>
   * - seuilcond
     - Threshold used for cond matrix
     - <class 'float'>
   * - NITT
     - Number of iterations
     - <class 'int'>
   * - FITANGLE
     -
     - <class 'int'>
   * - FITPOLEFF
     -
     - <class 'int'>
   * - saveCOV
     -
     - <class 'int'>
   * - saveCO
     -
     - <class 'int'>
   * - delta_psi
     -  Angle in degrees to add to all bolometers psi pointing ({dbpath}/{pixname}_REP6_ptg_TUPLE_2) to simulate an entire focal plane rotation.Omit or set to 0.0 for default/legacy behavior (no focal plane rotation)
     - <class 'float'>
   * - verbose
     -
     - <class 'int'>
   * - TEMPLATE_NSIDE
     -
     - <class 'int'>
   * - BUILDTF
     -
     - <class 'int'>
   * - SEED
     -
     - <class 'list'>
   * - GAINSTEP
     - Number of gain values to fit over the whole mission.
       If set to 1, only one gain value will be fitted, thus no variable gains.
       RD12 values are GAINSTEP=128 for 100GHz to 217GHz and GAINSTEP=32 for 353GHz to 857GHz.
     - <class 'int'>
   * - NADU
     -
     - <class 'list'>
   * - NADUSTEP
     -
     - <class 'list'>
   * - DOCNN
     -
     - <class 'list'>
   * - val_mean
     -
     - <class 'list'>
   * - w_mean
     -
     - <class 'list'>
   * - do_mean
     -
     - <class 'list'>
   * - CALLCNN
     -
     - <class 'str'>
   * - CNN_START
     -
     - <class 'int'>
   * - CNN_XSIZE
     -
     - <class 'int'>
   * - CNN_YSIZE
     -
     - <class 'int'>
   * - CNN_RESIDU
     -
     - <class 'float'>
   * - CNN_WEIGHTS
     -
     - <class 'str'>
   * - N_IN_ITT
     -
     - <class 'int'>
   * - INST_CNN
     -
     - <class 'str'>
   * - MAP_CNN
     -
     - <class 'str'>
   * - Calibration
     -
     - <class 'list'>
   * - NEP
     -
     - <class 'list'>
   * - CrossPol
     - Table of Polarization efficiency [%] for each detectors
     - <class 'list'>
   * - Monop
     -
     - <class 'list'>
   * - OUT_NOPOL
     -
     - <class 'list'>
   * - bolomask
     -
     - <class 'list'>
   * - beg_surv
     -
     - <class 'list'>
   * - end_surv
     -
     - <class 'list'>
   * - name_surv
     -
     - <class 'list'>
   * - MAPRINGS
     -
     - <class 'list'>
   * - Mask
     -
     - <class 'str'>
   * - projectionType
     -
     - <class 'str'>
   * - in_template_map
     -
     - <class 'list'>
   * - Signal_noPS
     -
     - <class 'list'>
   * - ADU
     -
     - <class 'list'>
   * - Badring
     -
     - <class 'list'>
   * - HPR_Calib
     - List of HPRs containing the total dipole component that will be subtracted,not fitted, from Signal_noPS.
     - <class 'list'>
   * - Hit_noPS
     -
     - <class 'list'>
   * - Ptg_noPS
     -
     - <class 'list'>
   * - Theo_HPR
     -
     - <class 'list'>
   * - Theo_MAP
     -
     - <class 'list'>
   * - phase
     - Input phase of the HPR to be projected to maps. Used for advance denoising
     - <class 'list'>
   * - rgcnn
     -
     - <class 'list'>
   * - bolo_map
     -
     - <class 'list'>
   * - Out_MAP
     -
     - <class 'list'>
   * - Out_VEC
     -
     - <class 'list'>
   * - Nside
     -
     - <class 'int'>
