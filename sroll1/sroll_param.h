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
 
#ifndef _SROLL_PARAM_H_
#define _SROLL_PARAM_H_


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
   *
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
  
  /* ---- List SEED ----------------------------------------------------
   *
   */
  PIOINT *SEED;
  PIOLONG n_SEED; /* ==number_of_SEED */
  
  /* ---- Parameter FITTHETA -------------------------------------------
   *
   */
  PIOINT FITTHETA;
  
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
  
  /* ---- Parameter CALCODUST ------------------------------------------
   *
   */
  PIOINT CALCODUST;
  
  /* ---- Parameter DOMAXVRAIE -----------------------------------------
   *
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
  
  /* ---- Parameter TESTPOL --------------------------------------------
   *
   */
  PIOINT TESTPOL;
  
  /* ---- Parameter RSTEP ----------------------------------------------
   *
   */
  PIOINT RSTEP;
  
  /* ---- Parameter GAINSTEP -------------------------------------------
   *
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
  
  /* ---- Parameter AVGCO ----------------------------------------------
   *
   */
  PIODOUBLE AVGCO;
  
  /* ---- Parameter AVGDUST --------------------------------------------
   *
   */
  PIODOUBLE AVGDUST;
  
  /* ---- Parameter NADU -----------------------------------------------
   *
   */
  PIOINT NADU;
  
  /* ---- Parameter REMDIP ---------------------------------------------
   *
   */
  PIOINT REMDIP;

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
   * Input signal HPR in Watts
   * unused if stim_paramfiles is given (input HPR will then be taken
   * from stim ouptut)
   */
  PIOSTRING *Signal_noPS;
  PIOBYTE flag_Signal_noPS;  /* ==_PAR_TRUE if Signal_noPS is present */
  PIOLONG n_Signal_noPS; /* ==number_of_Signal_noPS */
  
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
  
  /* ---- Optional Input Theo_CO_Q -------------------------------------
   * Input CO Map (theoritical) refitted inside the software
   */
  PIOSTRING Theo_CO_Q;
  PIOBYTE flag_Theo_CO_Q;  /* ==_PAR_TRUE if Theo_CO_Q is present */
  
  /* ---- Optional Input Theo_CO_U -------------------------------------
   * Input CO Map (theoritical) refitted inside the software
   */
  PIOSTRING Theo_CO_U;
  PIOBYTE flag_Theo_CO_U;  /* ==_PAR_TRUE if Theo_CO_U is present */
  
  /* ---- Optional Input Theo_FREEFREE ---------------------------------
   * Input FREEFREE Map (theoritical) refitted inside the
   * software
   */
  PIOSTRING Theo_FREEFREE;
  PIOBYTE flag_Theo_FREEFREE;  /* ==_PAR_TRUE if Theo_FREEFREE is present */
  
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
   *
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
   * list of ring cuts (eg. surveys/hm/oddeven) used to produce Maps for each bolomask
   */
  PIOINT *MAPRINGS;
  PIOBYTE flag_MAPRINGS;  /* ==_PAR_TRUE if MAPRINGS is present */
  PIOLONG n_MAPRINGS; /* ==number_of_MAP */
  
  /* ---- Optional Parameter ADDDIP ------------------------------------
   * set to 1 to add the dipole HPR to the input HPR at the beginning of sroll
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
 
} sroll_parContent;

#endif
