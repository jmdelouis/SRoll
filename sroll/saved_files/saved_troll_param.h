/*
 *
 * Template for hprcalpol6 parameter file reader include file
 * Automatically produced from ./../xml/hprcalpol6.pie
 *
 * Created by delouis on Thu Mar 26 11:12:14 2015 (PIE version v00-05-67)
 *
 * Command line was something like :
 * PIE ./../xml/hprcalpol6.pie MYFILE.h -Lh
 *
 */

/*
 *
 * This include file contains the definition of a structure that will
 * be filled by reading the parameter file.
 *
 */

#ifndef _TROLL_PARAM_H_
#define _TROLL_PARAM_H_


#include "no_dmc_piolib_type_def.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


/*
 * parcontent : structure containing the read parameter file
 */

typedef struct {

  /* ---- Parameter BeginRing ------------------------------------------
   * First ring to consider (included)
   */
  PIOLONG BeginRing;

  /* ---- Parameter EndRing --------------------------------------------
   * Last ring to consider (included)
   */
  PIOLONG EndRing;

  /* ---- List CrossPol ------------------------------------------------
   *
   */
  PIODOUBLE *CrossPol;
  PIOLONG n_CrossPol; /* ==number_of_CrossPol */

  /* ---- List ADDCO13 ---------------------------------------------
   *
   */
  PIODOUBLE *ADDCO13;
  PIOBYTE flag_ADDCO13;  /* ==_PAR_TRUE if ADDCO13 is present */
  PIOLONG n_ADDCO13; /* ==number_of_ADDCO13 */

  /* ---- List Calibration ---------------------------------------------
   *
   */
  PIODOUBLE *Calibration;
  PIOLONG n_Calibration; /* ==number_of_Calibration */

  /* ---- Parameter seuilcond ------------------------------------------
   * Threshold used for cond matrix
   */
  PIODOUBLE seuilcond;

  /* ---- List NEP -----------------------------------------------------
   * http://wiki.planck.fr/index.php/Proc/DetectorWeightingInMaps
   */
  PIODOUBLE *NEP;
  PIOLONG n_NEP; /* ==number_of_NEP */

  /* ---- List Monop ---------------------------------------------------
   *
   */
  PIODOUBLE *Monop;
  PIOLONG n_Monop; /* ==number_of_Monop */

  /* ---- List bolomask ------------------------------------------------
   *
   */
  PIOINT *bolomask;
  PIOLONG n_bolomask; /* ==number_of_bolomask */

  /* ---- List CNN ----------------------------------------------------
   *
   */
  PIOINT *DOCNN;
  PIOBYTE flag_DOCNN;  /* ==_PAR_TRUE if DOCNN is present */
  PIOLONG n_DOCNN; /* ==number_of_DOCNN */

  /* ---- List SEED ----------------------------------------------------
   *
   */
  PIOINT *SEED;
  PIOLONG n_SEED; /* ==number_of_SEED */

  /* ---- Parameter FITANGLE -------------------------------------------
   *
   */
  PIOINT FITANGLE;

  /* ---- Parameter FITPOLEFF -------------------------------------------
   *
   */
  PIOINT FITPOLEFF;

  /* ---- List OUT_NOPOL -----------------------------------------------
   *
   */
  PIOINT *OUT_NOPOL;
  PIOLONG n_OUT_NOPOL; /* ==number_of_OUT_NOPOL */

  /* ---- List FSLCOEF -------------------------------------------------
   *
   */
  PIOFLOAT *FSLCOEF;
  PIOLONG n_FSLCOEF; /* ==number_of_FSLCOEF */

  /* ---- Parameter D_NOPOL --------------------------------------------
   *
   */
  PIOINT D_NOPOL;

  /* ---- Parameter SAVEINTMAP -----------------------------------------
   *
   */
  PIOINT SAVEINTMAP;

  /* ---- Parameter CALLCNN -----------------------------------------
   *
   */
  PIOSTRING CALLCNN  ;
  PIOBYTE flag_CALLCNN;  /* ==_PAR_TRUE if CALLCNN is present */

  /* ---- Parameter CNN_ITT -----------------------------------------
   *
   */
  PIOINT CNN_ITT;
  PIOBYTE flag_CNN_ITT;  /* ==_PAR_TRUE if CNN_ITT is present */

  /* ---- Parameter CNN_ITERID -----------------------------------------
   *
   */
  PIOINT CNN_ITERID;
  PIOBYTE flag_CNN_ITERID;  /* ==_PAR_TRUE if CNN_ITERID is present */

  /* ---- Parameter CNN_TMPID -----------------------------------------
   *
   */
  PIOSTRING CNN_TMPID;
  PIOBYTE flag_CNN_TMPID;  /* ==_PAR_TRUE if CNN_TMPID is present */

  /* ---- Parameter PHASECNN -----------------------------------------
   *
   */
  PIOINT PHASECNN;
  PIOBYTE flag_PHASECNN;  /* ==_PAR_TRUE if PHASECNN is present */

  /* ---- Parameter CNN_RESIDU -----------------------------------------
   *
   */
  PIOFLOAT CNN_RESIDU;
  PIOBYTE flag_CNN_RESIDU;  /* ==_PAR_TRUE if CNN_RESIDU is present */

  /* ---- Parameter CNN_XSIZE -----------------------------------------
   *
   */
  PIOINT CNN_XSIZE;
  PIOBYTE flag_CNN_XSIZE;  /* ==_PAR_TRUE if CNN_XSIZE is present */

  /* ---- Parameter CNN_YSIZE -----------------------------------------
   *
   */
  PIOINT CNN_YSIZE;
  PIOBYTE flag_CNN_YSIZE;  /* ==_PAR_TRUE if CNN_YSIZE is present */

  /* ---- Parameter CNN_WEIGHTS -----------------------------------------
   *
   */
  PIOSTRING CNN_WEIGHTS;
  PIOBYTE flag_CNN_WEIGHTS;  /* ==_PAR_TRUE if CNN_WEIGHTS is present */

  /* ---- Parameter BUILDTF -----------------------------------------
   *
   */
  PIOINT  BUILDTF;
  PIOBYTE flag_BUILDTF;  /* ==_PAR_TRUE if BUILDTF is present */

  /* ---- Parameter AVFDUST -----------------------------------------
   *
   */
  PIOINT AVFDUST;

  /* ---- Parameter AVFCO -----------------------------------------
   *
   */
  PIOINT AVFCO;

  /* ---- Parameter AVFSYNC -----------------------------------------
   *
   */
  PIOINT AVFSYNC;

  /* ---- Parameter AVFFF -----------------------------------------
   *
   */
  PIOINT AVFFF;

  /* ---- Parameter CALCODUST ------------------------------------------
   *
   */
  PIOINT CALCODUST;

  /* ---- Parameter DOMAXVRAIE -----------------------------------------
   * if set to 0, maps for all surveys and half missions and odd and even rings
   * will be produced, for each map prefix in the MAP parameter list.
   * If set to 1, only one map from BeginRing to EndRing will be produced, per map prefix.
   */
  PIOINT DOMAXVRAIE;

  /* ---- Parameter DOGAINDIP ------------------------------------------
   *
   */
  PIOINT DOGAINDIP;

  /* ---- Parameter CUTRG ----------------------------------------------
   *
   */
  PIOINT CUTRG;

  /* ---- Optional Parameter AVGR0 -------------------------------------
   *
   */
  PIODOUBLE AVGR0;
  PIOBYTE flag_AVGR0;  /* ==_PAR_TRUE if AVGR0 is present */

  /* ---- Parameter DODIPCAL -------------------------------------------
   *
   */
  PIOINT DODIPCAL;

  /* ---- Optional Parameter TEMPLATE_NSIDE -------------------------------------------
   *
   */
  PIOINT TEMPLATE_NSIDE;
  PIOBYTE flag_TEMPLATE_NSIDE;  /* ==_PAR_TRUE if TEMPLATE_NSIDE is present */
  
  /* ---- Parameter TESTPOL --------------------------------------------
   * if set to 0, activates the sroll simulation code (without stim),
   * where colored noise and ADCNL residuals are added to the input signal HPR.
   * If set to 4, only the ADCNL residuals part is added (?). Other values allowed
   * in the code are: 3 (commented), 7 (?). Ignored (forced to -1) if stim_paramfiles is set.
   * ** SHOULD BE SET TO 1 **
   */
  PIOINT TESTPOL;

  /* ---- Parameter RSTEP ----------------------------------------------
   * ring step. If set to 1, all the input rings are processed by sroll,
   * if set to 10, 1 ring out of 10 is used, and so on. Useful to quick check
   * some parameter combinations, as run time is significantly reduced.
   * Use RSTEP=1 for actual productions.
   */
  PIOINT RSTEP;

  /* ---- Parameter GAINSTEP -------------------------------------------
   * number of gain values to fit over the whole mission.
   * If set to 1, only one gain value will be fitted, thus no variable gains.
   * RD12 values are GAINSTEP=128 for 100GHz to 217GHz and GAINSTEP=32 for 353GHz to 857GHz.
   */
  PIOINT GAINSTEP;

  /* ---- Parameter Nside ----------------------------------------------
   *
   */
  PIOINT Nside;

  /* ---- Parameter XI2STOP --------------------------------------------
   *
   */
  PIODOUBLE XI2STOP;


  /* ---- Parameter AVGDUST --------------------------------------------
   *
   */
  PIODOUBLE AVGDUST;

  /* ---- Parameter NITT -----------------------------------------------
   *
   */
  PIOINT NITT;

  /* ---- Parameter NADU -----------------------------------------------
   *
   */
  PIOINT *NADU;
  PIOINT n_NADU; /* ==number_NADU */

  /* ---- Parameter NADUSTEP -----------------------------------------------
   *
   */
  PIOINT *NADUSTEP;
  PIOINT n_NADUSTEP; /* ==number_NADUSTEP */

  /* ---- Parameter REMDIP ---------------------------------------------
   *
   */
  PIOINT REMDIP;

  /* ---- Parameter in_template_map_I ---------------------------------------------
   *
   */
  PIOSTRING *in_template_map_I;
  PIOBYTE flag_in_template_map_I;  /* ==_PAR_TRUE if in_template_map_I is present */
  PIOLONG n_in_template_map_I; /* ==number_in_template_map_I */

  /* ---- Parameter in_template_map_Q ---------------------------------------------
   *
   */
  PIOSTRING *in_template_map_Q;
  PIOBYTE flag_in_template_map_Q;  /* ==_PAR_TRUE if in_template_map_Q is present */
  PIOLONG n_in_template_map_Q; /* ==number_in_template_map_Q */

  /* ---- Parameter in_template_map_U ---------------------------------------------
   *
   */
  PIOSTRING *in_template_map_U;
  PIOBYTE flag_in_template_map_U;  /* ==_PAR_TRUE if in_template_map_U is present */
  PIOLONG n_in_template_map_U; /* ==number_in_template_map_U */


  /* ---- Optional Parameter REFMAPI --------------------------------------
   *
   */
  PIOSTRING REFMAPI;
  PIOBYTE flag_REFMAPI;  /* ==_PAR_TRUE if REFMAPI is present */

  /* ---- Optional Parameter REFMAPQ --------------------------------------
   *
   */
  PIOSTRING REFMAPQ;
  PIOBYTE flag_REFMAPQ;  /* ==_PAR_TRUE if REFMAPQ is present */

  /* ---- Optional Parameter REFMAPU --------------------------------------
   *
   */
  PIOSTRING REFMAPU;
  PIOBYTE flag_REFMAPU;  /* ==_PAR_TRUE if REFMAPU is present */

  /* ---- Optional Parameter REMHDIP ---------------------------------------------
   *
   */
  PIOINT REMHDIP;
  PIOBYTE flag_REMHDIP;  /* ==_PAR_TRUE if REMHDIP is present */

  /* ---- Optional Input TEMPLATEMAP ---------------------------------------
   *
   */
  PIOSTRING TEMPLATEMAP;
  PIOBYTE flag_TEMPLATEMAP;  /* ==_PAR_TRUE if TEMPLATEMAP is present */

  /* ---- Optional InputList Signal_noPS -------------------------------
   * Input signal HPR to be projected to maps.
   * Signal_noPS units multiplied by Calibration must be KCMB.
   * unused if stim_paramfiles is given (input HPR will then be taken from stim ouptut)
   */
  PIOSTRING *Signal_noPS;
  PIOBYTE flag_Signal_noPS;  /* ==_PAR_TRUE if Signal_noPS is present */
  PIOLONG n_Signal_noPS; /* ==number_of_Signal_noPS */

  /* ---- Optional InputList phase -------------------------------
   * Input phase of the HPR to be projected to maps. Used for advance denoising
   */
  PIOSTRING *phase;
  PIOBYTE flag_phase;  /* ==_PAR_TRUE if phase is present */
  PIOLONG n_phase; /* ==number_of_phase */

  /* ---- Optional InputList rgcnn -------------------------------
   * Input phase of the HPR to be projected to maps. Used for advance denoising
   */
  PIOSTRING *rgcnn;
  PIOBYTE flag_rgcnn;  /* ==_PAR_TRUE if rgcnn is present */
  PIOLONG n_rgcnn; /* ==number_of_rgcnn */

  /* ---- Optional InputList invgi -------------------------------------
   *
   */
  PIOSTRING *invgi;
  PIOBYTE flag_invgi;  /* ==_PAR_TRUE if invgi is present */
  PIOLONG n_invgi; /* ==number_of_invgi */

  /* ---- Optional InputList fsl ---------------------------------------
   *
   */
  PIOSTRING *fsl;
  PIOBYTE flag_fsl;  /* ==_PAR_TRUE if fsl is present */
  PIOLONG n_fsl; /* ==number_of_fsl */

  /* ---- Optional InputList Theo_noPS ---------------------------------
   *
   */
  PIOSTRING *Theo_noPS;
  PIOBYTE flag_Theo_noPS;  /* ==_PAR_TRUE if Theo_noPS is present */
  PIOLONG n_Theo_noPS; /* ==number_of_Theo_noPS */

  /* ---- Optional InputList ADU ---------------------------------------
   *
   */
  PIOSTRING *ADU;
  PIOBYTE flag_ADU;  /* ==_PAR_TRUE if ADU is present */
  PIOLONG n_ADU; /* ==number_of_ADU */

  /* ---- Optional Input Theo_CO ---------------------------------------
   * Input CO Map (theoritical) refitted inside the software
   */
  PIOSTRING Theo_CO;
  PIOBYTE flag_Theo_CO;  /* ==_PAR_TRUE if Theo_CO is present */

  /* ---- Optional Input Theo_13CO ---------------------------------------
   * Input 13CO Map (theoritical) refitted inside the software
   */
  PIOSTRING Theo_13CO;
  PIOBYTE flag_Theo_13CO;  /* ==_PAR_TRUE if Theo_13CO is present */

  /* ---- Optional Input in_polar_fit_Q -------------------------------------
   * Input CO Map (theoritical) refitted inside the software
   */
  PIOSTRING *in_polar_fit_Q;
  PIOBYTE flag_in_polar_fit_Q;  /* ==_PAR_TRUE if in_polar_fit_Q is present */
  PIOLONG n_in_polar_fit_Q; /* ==number_in_polar_fit_Q */

  /* ---- Optional Input in_polar_fit_U -------------------------------------
   * Input CO Map (theoritical) refitted inside the software
   */
  PIOSTRING *in_polar_fit_U;
  PIOBYTE flag_in_polar_fit_U;  /* ==_PAR_TRUE if in_polar_fit_U is present */
  PIOLONG n_in_polar_fit_U; /* ==number_in_polar_fit_U */

  /* ---- Optional Input Theo_FREEFREE ---------------------------------
   * Input FREEFREE Map (theoritical) refitted inside the
   * software
   */
  PIOSTRING Theo_FREEFREE;
  PIOBYTE flag_Theo_FREEFREE;  /* ==_PAR_TRUE if Theo_FREEFREE is present */

  /* ---- Optional Input in_synchro_map_I -----------------------------------
   * Input Dust Map (theoritical) refitted inside the software
   */
  PIOSTRING in_synchro_map_I;
  PIOBYTE flag_in_synchro_map_I;  /* ==_PAR_TRUE if in_synchro_map_I is present */

  /* ---- Optional Input in_synchro_map_Q -----------------------------------
   *  Input Dust Map (theoritical) refitted inside the software
   */
  PIOSTRING in_synchro_map_Q;
  PIOBYTE flag_in_synchro_map_Q;  /* ==_PAR_TRUE if in_synchro_map_Q is present */

  /* ---- Optional Input in_synchro_map_U -----------------------------------
   *  Input Dust Map (theoritical) refitted inside the software
   */
  PIOSTRING in_synchro_map_U;
  PIOBYTE flag_in_synchro_map_U;  /* ==_PAR_TRUE if in_synchro_map_U is present */

  /* ---- Optional Input Theo_Dust_I -----------------------------------
   * Input Dust Map (theoritical) refitted inside the software
   */
  PIOSTRING Theo_Dust_I;
  PIOBYTE flag_Theo_Dust_I;  /* ==_PAR_TRUE if Theo_Dust_I is present */

  /* ---- Optional Input Theo_Dust_Q -----------------------------------
   *  Input Dust Map (theoritical) refitted inside the software
   */
  PIOSTRING Theo_Dust_Q;
  PIOBYTE flag_Theo_Dust_Q;  /* ==_PAR_TRUE if Theo_Dust_Q is present */

  /* ---- Optional Input Theo_Dust_U -----------------------------------
   *  Input Dust Map (theoritical) refitted inside the software
   */
  PIOSTRING Theo_Dust_U;
  PIOBYTE flag_Theo_Dust_U;  /* ==_PAR_TRUE if Theo_Dust_U is present */

  /* ---- Optional Input Theo_Tdust_I -----------------------------------
   * Input Tdust Map (theoritical) refitted inside the software
   */
  PIOSTRING Theo_TDust_I;
  PIOBYTE flag_Theo_TDust_I;  /* ==_PAR_TRUE if Theo_Tdust_I is present */

  /* ---- Optional Input Theo_Tdust_Q -----------------------------------
   *  Input Tdust Map (theoritical) refitted inside the software
   */
  PIOSTRING Theo_TDust_Q;
  PIOBYTE flag_Theo_TDust_Q;  /* ==_PAR_TRUE if Theo_Tdust_Q is present */

  /* ---- Optional Input Theo_Tdust_U -----------------------------------
   *  Input Tdust Map (theoritical) refitted inside the software
   */
  PIOSTRING Theo_TDust_U;
  PIOBYTE flag_Theo_TDust_U;  /* ==_PAR_TRUE if Theo_Tdust_U is present */

  /* ---- Optional Input Sim_Dust_Q ------------------------------------
   *
   */
  PIOSTRING Sim_Dust_Q;
  PIOBYTE flag_Sim_Dust_Q;  /* ==_PAR_TRUE if Sim_Dust_Q is present */

  /* ---- Optional Input Sim_Dust_U ------------------------------------
   *
   */
  PIOSTRING Sim_Dust_U;
  PIOBYTE flag_Sim_Dust_U;  /* ==_PAR_TRUE if Sim_Dust_U is present */

  /* ---- Optional InputList Ptg_noPS ----------------------------------
   *
   */
  PIOSTRING *Ptg_noPS;
  PIOBYTE flag_Ptg_noPS;  /* ==_PAR_TRUE if Ptg_noPS is present */
  PIOLONG n_Ptg_noPS; /* ==number_of_Ptg_noPS */

  /* ---- Optional InputList Hit_noPS ----------------------------------
   *
   */
  PIOSTRING *Hit_noPS;
  PIOBYTE flag_Hit_noPS;  /* ==_PAR_TRUE if Hit_noPS is present */
  PIOLONG n_Hit_noPS; /* ==number_of_Hit_noPS */

  /* ---- Optional InputList ALMMAP ------------------------------------
   *
   */
  PIOSTRING *ALMMAP;
  PIOBYTE flag_ALMMAP;  /* ==_PAR_TRUE if ALMMAP is present */
  PIOLONG n_ALMMAP; /* ==number_of_ALMMAP */

  /* ---- Optional InputList DipOrb_noPS -------------------------------
   * list of HPRs containing the total dipole component that will be subtracted,
   * not fitted, from Signal_noPS. Must be KCMB.
   */
  PIOSTRING *DipOrb_noPS;
  PIOBYTE flag_DipOrb_noPS;  /* ==_PAR_TRUE if DipOrb_noPS is present */
  PIOLONG n_DipOrb_noPS; /* ==number_of_DipOrb_noPS */

  /* ---- Optional InputList Badring -----------------------------------
   *
   */
  PIOSTRING *Badring;
  PIOBYTE flag_Badring;  /* ==_PAR_TRUE if Badring is present */
  PIOLONG n_Badring; /* ==number_of_Badring */

  /* ---- Optional InputList VarGain -----------------------------------
   *
   */
  PIOSTRING *VarGain;
  PIOBYTE flag_VarGain;  /* ==_PAR_TRUE if VarGain is present */
  PIOLONG n_VarGain; /* ==number_of_VarGain */

  /* ---- Optional Input Mask ------------------------------------------
   *
   */
  PIOSTRING Mask;
  PIOBYTE flag_Mask;  /* ==_PAR_TRUE if Mask is present */

  /* ---- Optional OutputList MAP --------------------------------------
   *
   */
  PIOSTRING *MAP;
  PIOBYTE flag_MAP;  /* ==_PAR_TRUE if MAP is present */
  PIOLONG n_MAP; /* ==number_of_MAP */

  /* ---- Optional OutputList Out_Offset -------------------------------
   *
   */
  PIOSTRING *Out_Offset;
  PIOBYTE flag_Out_Offset;  /* ==_PAR_TRUE if Out_Offset is present */
  PIOLONG n_Out_Offset; /* ==number_of_Out_Offset */

  /* ---- Optional OutputList Out_Offset_corr --------------------------
   *
   */
  PIOSTRING *Out_Offset_corr;
  PIOBYTE flag_Out_Offset_corr;  /* ==_PAR_TRUE if Out_Offset_corr is present */
  PIOLONG n_Out_Offset_corr; /* ==number_of_Out_Offset_corr */

  /* ---- Optional OutputList Out_xi2 ----------------------------------
   *
   */
  PIOSTRING *Out_xi2;
  PIOBYTE flag_Out_xi2;  /* ==_PAR_TRUE if Out_xi2 is present */
  PIOLONG n_Out_xi2; /* ==number_of_Out_xi2 */

  /* ---- Optional OutputList Out_xi2_corr -----------------------------
   *
   */
  PIOSTRING *Out_xi2_corr;
  PIOBYTE flag_Out_xi2_corr;  /* ==_PAR_TRUE if Out_xi2_corr is present */
  PIOLONG n_Out_xi2_corr; /* ==number_of_Out_xi2_corr */

  /* ---- Optional Parameter verbose -----------------------------------
   * Verbosity level. 0 : normal, 1:verbose
   */
  PIOINT verbose;
  PIOBYTE flag_verbose;  /* ==_PAR_TRUE if verbose is present */

  /* ---- Optional Parameter dmc_output_path ---------------------------
   * Special slot to store the output logging path
   */
  PIOSTRING dmc_output_path;
  PIOBYTE flag_dmc_output_path;  /* ==_PAR_TRUE if dmc_output_path is present */

  /* ---- Optional Parameter dmc_error_path ----------------------------
   * Special slot to store the error logging path
   */
  PIOSTRING dmc_error_path;
  PIOBYTE flag_dmc_error_path;  /* ==_PAR_TRUE if dmc_error_path is present */



  /* ---- Optional InputList stim_paramfiles -------------------------------
   * Parameter files for TOI simulations, one per bolo, used if TESTPOL==0
   */
  PIOSTRING *stim_paramfiles;
  PIOBYTE flag_stim_paramfiles;  /* ==_PAR_TRUE if stim_paramfiles is present */
  PIOLONG n_stim_paramfiles; /* ==number_of_stim_paramfiles */

  /* ---- Optional OutputList MAPRINGS --------------------------------------
   * Maps to produce for each bolomask
   */
  PIOINT *MAPRINGS;
  PIOBYTE flag_MAPRINGS;  /* ==_PAR_TRUE if MAPRINGS is present */
  PIOLONG n_MAPRINGS; /* ==number_of_MAP */

  /* ---- Optional Parameter ADDDIP ------------------------------------
   * set to 1 to add the dipole HPR to the input HPR at the beginning of troll
   */
  PIOINT ADDDIP;
  PIOBYTE flag_ADDDIP;

  /* ---- Optional Parameter KCMBIN ------------------------------------
   * set to 1 if Signal_noPS input HPR is in KCMB instead of Watts
   * when KCMBIN=1, Signal_noPs is first converted to Watts using the Calibration parameter
   */
  PIOINT KCMBIN;
  PIOBYTE flag_KCMBIN;

  /* ---- Optional InputList addHPR_name -------------------------------
   * list of HPR object names of type PIOFLOAT to add to Signal_noPS input HPRs
   * must be a multiple of number of bolometers
   */
  PIOSTRING *addHPR_name;
  PIOBYTE   flag_addHPR_name;
  PIOLONG   n_addHPR_name;

  /* ---- Optional InputList addHPR_factor -------------------------------
   * list of multiplicative factor for each addHPR_name
   * must have the same number of elements as addHPR_name
   * if only one value is given, it will be used for all addHPR_name objects
   * default is 1.0 (simply add the addHPR_name to Signal_noPS)
   */
  PIOFLOAT *addHPR_factor;
  PIOBYTE  flag_addHPR_factor;
  PIOLONG  n_addHPR_factor;

  /* ---- Optional InputList addHPR_watts -------------------------------
   * if set to 1, the corresponding addHPR object is in Watts and will
   * be converted to KCMB before being used
   * if only one value is given, it will be used for all addHPR_name objects
   * default is 0 (HPRs in KCMB)
   */
  PIOINT  *addHPR_watts;
  PIOBYTE flag_addHPR_watts;
  PIOLONG n_addHPR_watts;

  /* ---- Optional Parameter AVG12CO ----------------------------------
   * when using bolometers from different frequencies,
   * if AVG12CO is given, 12CO contribution will be forced to AVG12CO value
   * for 143GHz bolometer and free for the others
   * eg: AVG12CO = 0.0
   */
  PIODOUBLE AVG12CO;
  PIOBYTE flag_AVG12CO;  /* = 1 if AVG12CO is present */
  PIODOUBLE AVG13CO;
  PIOBYTE flag_AVG13CO;  /* = 1 if AVG12CO is present */

  /* ---- Optional Parameter AVGDUST100 --------------------------------
   * when using bolometers from different frequencies,
   * if AVGDUST100 is given, DUST contribution will be forced to AVGDUST100 value
   * for 100GHz bolometer and free for the others
   * eg: AVGDUST100 = 0.02
   */
  PIODOUBLE AVGDUST100;
  PIOBYTE flag_AVGDUST100;  /* = 1 if AVGDUST100 is present */

  PIODOUBLE AVGFREEFREE;
  PIOBYTE flag_AVGFREEFREE;
  PIODOUBLE AVGSYNCHRO;
  PIOBYTE flag_AVGSYNCHRO;

  /* ---- Optional Parameter saveCOV ---------------------------------------------
   * if set to 0, inverse covariance matrix (II, IQ, IU, QQ, QU, UU) and HIT count won't be saved
   * if omitted or set to 1 (default), inverse cov and hit will be saved
   */
  PIOINT saveCOV;
  PIOBYTE flag_saveCOV;  /* ==_PAR_TRUE if saveCOV is present */

  /* ---- Optional Parameter saveCO ---------------------------------------------
   * if set to 1, monobolometer CO maps will be saved
   * if omitted or set to 0 (default), no monobolo CO maps will be saved
   */
  PIOINT saveCO;
  PIOBYTE flag_saveCO;  /* ==_PAR_TRUE if saveCO is present */

  /* ---- Optional Parameter delta_psi ---------------------------------------------
   * angle in degrees to add to all bolometers psi pointing ({dbpath}/{pixname}_REP6_ptg_TUPLE_2)
   * to simulate an entire focal plane rotation.
   * omit or set to 0.0 for default/legacy behavior (no focal plane rotation)
   */
  PIOFLOAT delta_psi;
  PIOBYTE flag_delta_psi;  /* = 1 if delta_psi is present */

} troll_parContent;

#endif
