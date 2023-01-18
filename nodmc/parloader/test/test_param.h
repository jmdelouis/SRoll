/*
 * 
 * This include file contains the definition of a structure that will
 * be filled by reading a test parameter file.
 *
 */
 
#ifndef _TEST_PARAM_H_
#define _TEST_PARAM_H_


#include "no_dmc_piolib_type_def.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


/*
 * test_parContent : structure containing the read parameter file
 */
 
typedef struct {
  
  /* ---- Parameter BeginRing ------------------------------------------
   *  first ring to consider (included)
   */
  PIOLONG BeginRing;
  
  /* ---- Parameter EndRing --------------------------------------------
   *  last ring to consider (included)
   */
  PIOLONG EndRing;
  
  /* ---- List CrossPol ------------------------------------------------
   *  last ring to consider (included)
   */
  PIODOUBLE *CrossPol;
  PIOLONG n_CrossPol; /* ==number_of_CrossPol */
  
  /* ---- List Calibration ---------------------------------------------
   *  last ring to consider (included)
   */
  PIODOUBLE *Calibration;
  PIOLONG n_Calibration; /* ==number_of_Calibration */

  /* ---- List NEP -----------------------------------------------------
   *  last ring to consider (included)
   */
  PIODOUBLE *NEP;
  PIOLONG n_NEP; /* ==number_of_NEP */
  
  /* ---- List Monop ---------------------------------------------------
   *  last ring to consider (included)
   */
  PIODOUBLE *Monop;
  PIOLONG n_Monop; /* ==number_of_Monop */
  
  /* ---- List bolomask ------------------------------------------------
   *  last ring to consider (included)
   */
  PIOINT *bolomask;
  PIOLONG n_bolomask; /* ==number_of_bolomask */
  
  /* ---- List SEED ----------------------------------------------------
   *  TESTPOL
   */
  PIOINT *SEED;
  PIOLONG n_SEED; /* ==number_of_SEED */
  
  /* ---- Parameter FITTHETA -------------------------------------------
   *  FITTHETA
   */
  PIOINT FITTHETA;
  
  /* ---- List OUT_NOPOL -----------------------------------------------
   *  FITTHETA
   */
  PIOINT *OUT_NOPOL;
  PIOLONG n_OUT_NOPOL; /* ==number_of_OUT_NOPOL */
  
  /* ---- List FSLCOEF -------------------------------------------------
   *  FITTHETA
   */
  PIOFLOAT *FSLCOEF;
  PIOLONG n_FSLCOEF; /* ==number_of_FSLCOEF */
  
  /* ---- Parameter D_NOPOL --------------------------------------------
   *  FITTHETA
   */
  PIOINT D_NOPOL;
  
  /* ---- Parameter SAVEINTMAP -----------------------------------------
   *  TESTPOL
   */
  PIOINT SAVEINTMAP;
  
  /* ---- Parameter CALCODUST ------------------------------------------
   *  TESTPOL
   */
  PIOINT CALCODUST;
  
  /* ---- Parameter DOMAXVRAIE -----------------------------------------
   *  TESTPOL
   */
  PIOINT DOMAXVRAIE;
  
  /* ---- Parameter DOGAINDIP ------------------------------------------
   *  TESTPOL
   */
  PIOINT DOGAINDIP;
  
  /* ---- Parameter CUTRG ----------------------------------------------
   *  TESTPOL
   */
  PIOINT CUTRG;
  
  /* ---- Optional Parameter AVGR0 -------------------------------------
   *  TESTPOL
   */
  PIODOUBLE AVGR0;
  PIOBYTE flag_AVGR0;  /* ==_PAR_TRUE if AVGR0 is present */
  
  /* ---- Parameter DODIPCAL -------------------------------------------
   *  TESTPOL
   */
  PIOINT DODIPCAL;
  
  /* ---- Parameter TESTPOL --------------------------------------------
   *  TESTPOL
   */
  PIOINT TESTPOL;
  
  /* ---- Parameter RSTEP ----------------------------------------------
   *  TESTPOL
   */
  PIOINT RSTEP;
  
  /* ---- Parameter GAINSTEP -------------------------------------------
   *  TESTPOL
   */
  PIOINT GAINSTEP;
  
  /* ---- Parameter Nside ----------------------------------------------
   *  Nside
   */
  PIOINT Nside;
  
  /* ---- Parameter XI2STOP --------------------------------------------
   *  XI2STOP
   */
  PIODOUBLE XI2STOP;
  
  /* ---- Parameter AVGCO ----------------------------------------------
   *  XI2STOP
   */
  PIODOUBLE AVGCO;
  
  /* ---- Parameter AVGDUST --------------------------------------------
   *  XI2STOP
   */
  PIODOUBLE AVGDUST;
  
  /* ---- Parameter DEGREE ---------------------------------------------
   *  XI2STOP
   */
  PIODOUBLE DEGREE;
  
  /* ---- Parameter NADU -----------------------------------------------
   *  XI2STOP
   */
  PIOINT NADU;
  
  /* ---- Parameter REMDIP ---------------------------------------------
   *  XI2STOP
   */
  PIOINT REMDIP;
  
  /* ---- Optional InputList Signal_noPS -------------------------------
   *  input timeline RAW timelines : if it is not there compute
   *  the dnl correction.
   */
  PIOSTRING *Signal_noPS;
  PIOBYTE flag_Signal_noPS;  /* ==_PAR_TRUE if Signal_noPS is present */
  PIOLONG n_Signal_noPS; /* ==number_of_Signal_noPS */
  
  /* ---- Optional InputList invgi -------------------------------------
   *  input timeline RAW timelines : if it is not there compute
   *  the dnl correction.
   */
  PIOSTRING *invgi;
  PIOBYTE flag_invgi;  /* ==_PAR_TRUE if invgi is present */
  PIOLONG n_invgi; /* ==number_of_invgi */
  
  /* ---- Optional InputList fsl ---------------------------------------
   *  input timeline RAW timelines : if it is not there compute
   *  the dnl correction.
   */
  PIOSTRING *fsl;
  PIOBYTE flag_fsl;  /* ==_PAR_TRUE if fsl is present */
  PIOLONG n_fsl; /* ==number_of_fsl */
  
  /* ---- Optional InputList Theo_noPS ---------------------------------
   *  input timeline RAW timelines : if it is not there compute
   *  the dnl correction.
   */
  PIOSTRING *Theo_noPS;
  PIOBYTE flag_Theo_noPS;  /* ==_PAR_TRUE if Theo_noPS is present */
  PIOLONG n_Theo_noPS; /* ==number_of_Theo_noPS */
  
  /* ---- Optional InputList ADU ---------------------------------------
   *  input timeline RAW timelines : if it is not there compute
   *  the dnl correction.
   */
  PIOSTRING *ADU;
  PIOBYTE flag_ADU;  /* ==_PAR_TRUE if ADU is present */
  PIOLONG n_ADU; /* ==number_of_ADU */
  
  /* ---- Optional Input Theo_CO ---------------------------------------
   *  Input CO Map (theoritical) refitted inside the software
   */
  PIOSTRING Theo_CO;
  PIOBYTE flag_Theo_CO;  /* ==_PAR_TRUE if Theo_CO is present */
  
  /* ---- Optional Input Theo_CO_Q -------------------------------------
   *  Input CO Map (theoritical) refitted inside the software
   */
  PIOSTRING Theo_CO_Q;
  PIOBYTE flag_Theo_CO_Q;  /* ==_PAR_TRUE if Theo_CO_Q is present */
  
  /* ---- Optional Input Theo_CO_U -------------------------------------
   *  Input CO Map (theoritical) refitted inside the software
   */
  PIOSTRING Theo_CO_U;
  PIOBYTE flag_Theo_CO_U;  /* ==_PAR_TRUE if Theo_CO_U is present */
  
  /* ---- Optional Input Theo_FREEFREE ---------------------------------
   *  Input FREEFREE Map (theoritical) refitted inside the
   *  software
   */
  PIOSTRING Theo_FREEFREE;
  PIOBYTE flag_Theo_FREEFREE;  /* ==_PAR_TRUE if Theo_FREEFREE is present */
  
  /* ---- Optional Input Theo_Dust_I -----------------------------------
   *  Input CO Map (theoritical) refitted inside the software
   */
  PIOSTRING Theo_Dust_I;
  PIOBYTE flag_Theo_Dust_I;  /* ==_PAR_TRUE if Theo_Dust_I is present */
  
  /* ---- Optional Input Theo_Dust_Q -----------------------------------
   *  Input CO Map (theoritical) refitted inside the software
   */
  PIOSTRING Theo_Dust_Q;
  PIOBYTE flag_Theo_Dust_Q;  /* ==_PAR_TRUE if Theo_Dust_Q is present */
  
  /* ---- Optional Input Theo_Dust_U -----------------------------------
   *  Input CO Map (theoritical) refitted inside the software
   */
  PIOSTRING Theo_Dust_U;
  PIOBYTE flag_Theo_Dust_U;  /* ==_PAR_TRUE if Theo_Dust_U is present */
  
  /* ---- Optional Input Sim_Dust_Q ------------------------------------
   *  Input CO Map (theoritical) refitted inside the software
   */
  PIOSTRING Sim_Dust_Q;
  PIOBYTE flag_Sim_Dust_Q;  /* ==_PAR_TRUE if Sim_Dust_Q is present */
  
  /* ---- Optional Input Sim_Dust_U ------------------------------------
   *  Input CO Map (theoritical) refitted inside the software
   */
  PIOSTRING Sim_Dust_U;
  PIOBYTE flag_Sim_Dust_U;  /* ==_PAR_TRUE if Sim_Dust_U is present */
  
  /* ---- Optional InputList Ptg_noPS ----------------------------------
   *  input timeline RAW timelines : if it is not there compute
   *  the dnl correction.
   */
  PIOSTRING *Ptg_noPS;
  PIOBYTE flag_Ptg_noPS;  /* ==_PAR_TRUE if Ptg_noPS is present */
  PIOLONG n_Ptg_noPS; /* ==number_of_Ptg_noPS */
  
  /* ---- Optional InputList Hit_noPS ----------------------------------
   *  input timeline RAW timelines : if it is not there compute
   *  the dnl correction.
   */
  PIOSTRING *Hit_noPS;
  PIOBYTE flag_Hit_noPS;  /* ==_PAR_TRUE if Hit_noPS is present */
  PIOLONG n_Hit_noPS; /* ==number_of_Hit_noPS */
  
  /* ---- Optional InputList ALMMAP ------------------------------------
   *  input timeline RAW timelines : if it is not there compute
   *  the dnl correction.
   */
  PIOSTRING *ALMMAP;
  PIOBYTE flag_ALMMAP;  /* ==_PAR_TRUE if ALMMAP is present */
  PIOLONG n_ALMMAP; /* ==number_of_ALMMAP */
  
  /* ---- Optional InputList DipOrb_noPS -------------------------------
   *  PARITY timeline
   */
  PIOSTRING *DipOrb_noPS;
  PIOBYTE flag_DipOrb_noPS;  /* ==_PAR_TRUE if DipOrb_noPS is present */
  PIOLONG n_DipOrb_noPS; /* ==number_of_DipOrb_noPS */
  
  /* ---- Optional InputList Badring -----------------------------------
   *  PARITY timeline
   */
  PIOSTRING *Badring;
  PIOBYTE flag_Badring;  /* ==_PAR_TRUE if Badring is present */
  PIOLONG n_Badring; /* ==number_of_Badring */
  
  /* ---- Optional InputList VarGain -----------------------------------
   *  PARITY timeline
   */
  PIOSTRING *VarGain;
  PIOBYTE flag_VarGain;  /* ==_PAR_TRUE if VarGain is present */
  PIOLONG n_VarGain; /* ==number_of_VarGain */
  
  /* ---- Optional Input Mask ------------------------------------------
   *  PARITY timeline
   */
  PIOSTRING Mask;
  PIOBYTE flag_Mask;  /* ==_PAR_TRUE if Mask is present */
  
  /* ---- Output MAPS1S2 -----------------------------------------------
   *  input timeline RAW timelines
   */
  PIOSTRING MAPS1S2;
  
  /* ---- Optional OutputList MAP --------------------------------------
   *  input timeline RAW timelines
   */
  PIOSTRING *MAP;
  PIOBYTE flag_MAP;  /* ==_PAR_TRUE if MAP is present */
  PIOLONG n_MAP; /* ==number_of_MAP */
  
  /* ---- Optional OutputList Out_Offset -------------------------------
   *  input timeline RAW timelines
   */
  PIOSTRING *Out_Offset;
  PIOBYTE flag_Out_Offset;  /* ==_PAR_TRUE if Out_Offset is present */
  PIOLONG n_Out_Offset; /* ==number_of_Out_Offset */
  
  /* ---- Optional OutputList Out_Offset_corr --------------------------
   *  input timeline RAW timelines
   */
  PIOSTRING *Out_Offset_corr;
  PIOBYTE flag_Out_Offset_corr;  /* ==_PAR_TRUE if Out_Offset_corr is present */
  PIOLONG n_Out_Offset_corr; /* ==number_of_Out_Offset_corr */
  
  /* ---- Optional OutputList Out_xi2 ----------------------------------
   *  input timeline RAW timelines
   */
  PIOSTRING *Out_xi2;
  PIOBYTE flag_Out_xi2;  /* ==_PAR_TRUE if Out_xi2 is present */
  PIOLONG n_Out_xi2; /* ==number_of_Out_xi2 */
  
  /* ---- Optional OutputList Out_xi2_corr -----------------------------
   *  input timeline RAW timelines
   */
  PIOSTRING *Out_xi2_corr;
  PIOBYTE flag_Out_xi2_corr;  /* ==_PAR_TRUE if Out_xi2_corr is present */
  PIOLONG n_Out_xi2_corr; /* ==number_of_Out_xi2_corr */
  
  /* ---- Optional Parameter verbose -----------------------------------
   *  Verbosity level. 0 : normal, 1:verbose
   */
  PIOINT verbose;
  PIOBYTE flag_verbose;  /* ==_PAR_TRUE if verbose is present */
  
  /* ---- Parameter mpi_Node -------------------------------------------
   *  Number of Node for mpi run
   */

  /* ---- Optional Parameter dmc_output_path ---------------------------
   *  Special slot to store the output logging path
   */
  PIOSTRING dmc_output_path;
  PIOBYTE flag_dmc_output_path;  /* ==_PAR_TRUE if dmc_output_path is present */
  
  /* ---- Optional Parameter dmc_error_path ----------------------------
   *  Special slot to store the error logging path
   */
  PIOSTRING dmc_error_path;
  PIOBYTE flag_dmc_error_path;  /* ==_PAR_TRUE if dmc_error_path is present */
  
  
 
} test_parContent;

#endif
