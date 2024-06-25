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

  /* ---- Parameter EndRing --------------------------------------------
   * Define the size of data to exchange between processor corresponding to a number of pixel 
   */
  PIOLONG MAXMPIBUFFER;

  
  /* ---- List ADDPOL ---------------------------------------------
   *
   */
  PIOSTRING ADDPOL;
  PIOBYTE flag_ADDPOL;  /* ==_PAR_TRUE if ADDPOL is present */
  
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

  /* ---- List CNN_CORE ---------------------------------------------
   *
   */
  PIOINT CNN_CORE;
  PIOBYTE flag_CNN_CORE;  /* ==_PAR_TRUE if CNN_CORE is present */
  
  /* ---- List SEED ----------------------------------------------------
   *
   */
  PIOINT *SEED;
  PIOLONG n_SEED; /* ==number_of_SEED */

  /* ---- List OUT_NOPOL -----------------------------------------------
   *
   */
  PIOINT *OUT_NOPOL;
  PIOLONG n_OUT_NOPOL; /* ==number_of_OUT_NOPOL */

  /* ---- List SUB_HPRCOEF -------------------------------------------------
   *
   */
  PIOFLOAT *SUB_HPRCOEF;
  PIOLONG n_SUB_HPRCOEF; /* ==number_of_SUB_HPRCOEF */

  /* ---- Parameter SAVEINTMAP -----------------------------------------
   *
   */
  PIOINT SAVEINTMAP;
  
  /* ---- Parameter CNN_LEARN_PARAM -----------------------------------------
   *
   */
  PIOINT CNN_LEARN_PARAM;
  PIOBYTE flag_CNN_LEARN_PARAM;  /* ==_PAR_TRUE if CNN_LEAN_PARAM is present */

  /* ---- Parameter CNN_TMPID -----------------------------------------
   *
   */
  PIOSTRING CNN_TMPID;
  PIOBYTE flag_CNN_TMPID;  /* ==_PAR_TRUE if CNN_TMPID is present */
 

  /*-----------------------------------------
   *
   */
  PIOINT  BUILDTF;
  PIOBYTE flag_BUILDTF;  /* ==_PAR_TRUE if BUILDTF is present */


  /* ---- Optional Parameter TEMPLATE_NSIDE -------------------------------------------
   *
   */
  PIOINT TEMPLATE_NSIDE;
  PIOBYTE flag_TEMPLATE_NSIDE;  /* ==_PAR_TRUE if TEMPLATE_NSIDE is present */

  /* ---- Parameter RSTEP ----------------------------------------------
   * ring step. If set to 1, all the input rings are processed by sroll,
   * if set to 10, 1 ring out of 10 is used, and so on. Useful to quick check
   * some parameter combinations, as run time is significantly reduced.
   * Use RSTEP=1 for actual productions.
   */
  PIOINT RSTEP;

  /* ---- Parameter RINGSIZE ----------------------------------------------
   * ring size in samples. 
   */
  PIOLONG RINGSIZE;

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
  /* ---- Parameter do_foscat ----------------------------------------------
   *
   */
  PIOINT * do_foscat;
  PIOLONG n_do_foscat;
  /* ---- Parameter do_templates ----------------------------------------------
   *
   */
  PIOINT * do_templates;
  PIOLONG n_do_templates;
  

  /* ---- Parameter NITT -----------------------------------------------
   *
   */
  PIOINT NITT;

  /* ---- Parameter N_IN_ITT -----------------------------------------------
   *
   */
  PIOINT N_IN_ITT;
  
  /* ---- Parameter S_IN_ITT -----------------------------------------------
   *
   */
  PIODOUBLE S_IN_ITT;
  
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

  /* ---- Parameter in_template_map ---------------------------------------------
   *
   */
  PIOSTRING *in_template_map;
  PIOBYTE flag_in_template_map;  /* ==_PAR_TRUE if in_template_map is present */
  PIOLONG n_in_template_map; /* ==number_in_template_map */

  /* ---- Parameter  number_val ---------------------------------------------
   *
   */  
  PIOINT number_val;
  
  
  /* ---- Parameter  beg_surv ---------------------------------------------
   *
   */  
  PIOLONG *beg_surv;
  PIOLONG n_beg_surv; /* ==number_of_beg_surv */
  
  /* ---- Parameter  end_surv ---------------------------------------------
   *
   */  
  PIOLONG *end_surv;
  PIOLONG n_end_surv; /* ==number_of_end_surv */
  
  /* ---- InputList name_surv ---------------------------------------------
   *
   */  
  PIOSTRING *name_surv;
  PIOLONG n_name_surv; /* ==number_of_name_surv */

  /* ---- Parameter  do_mean ---------------------------------------------
   *
   */  
  PIODOUBLE *do_mean;
  PIOLONG n_do_mean; /* ==number_of_do_mean */
 
  /* ---- Parameter  val_mean ---------------------------------------------
     *
     */  
    PIODOUBLE *val_mean;
    PIOLONG n_val_mean; /* ==number_of_do_mean */

  /* ---- Parameter  w_mean ---------------------------------------------
     *
     */  
    PIODOUBLE *w_mean;
    PIOLONG n_w_mean; /* ==number_of_do_mean */

   /* ---- Parameter  UNSSEN ---------------------------------------------
   *
   */  
    PIODOUBLE UNSEEN;

  /* ---- Optional Parameter projection ---------------------------------------------
   *
   */  
  PIOSTRING projection;
  PIOBYTE flag_projection;  /* ==_PAR_TRUE if projection is present */

  /* ---- Optional Parameter NORM_GAIN ---------------------------------------------
   *
   */
  PIOINT NORM_GAIN;
  
  /* ---- Optional Parameter REMOVE_CAL ---------------------------------------------
   *
   */
  PIOINT REMOVE_CAL;

  /* ---- Optional InputList Signal -------------------------------
   * Input signal HPR to be projected to maps.
   * Signal units multiplied by Calibration must be KCMB.
   * unused if stim_paramfiles is given (input HPR will then be taken from stim ouptut)
   */
  PIOSTRING *Signal;
  PIOBYTE flag_Signal;  /* ==_PAR_TRUE if Signal is present */
  PIOLONG n_Signal;     /* ==number_of_Signal */

  /* ---- Optional InputList External -------------------------------
   * Input External information of the HPR to be projected to maps. Used for advance denoising
   */
  PIOSTRING *External;
  PIOBYTE flag_External;  /* ==_PAR_TRUE if External is present */
  PIOLONG n_External;     /* ==number_of_External */

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
  PIOLONG n_invgi;     /* ==number_of_invgi */

  /* ---- Optional InputList Sub_HPR ---------------------------------------
   *
   */
  PIOSTRING *Sub_HPR;
  PIOBYTE flag_Sub_HPR;  /* ==_PAR_TRUE if Sub_HPR is present */
  PIOLONG n_Sub_HPR; /* ==number_of_Sub_HPR */
  
  /* ---- Optional Input SparseFunc ---------------------------------
   *
   */
  PIOSTRING SparseFunc;
  PIOBYTE flag_SparseFunc;  /* ==_PAR_TRUE if SparseFunc is present */

  /* ---- Optional Input DiagFunc ---------------------------------
   *
   */
  PIOSTRING DiagFunc;
  PIOBYTE flag_DiagFunc;  /* ==_PAR_TRUE if SparseFunc is present */

  /* ---- Optional InputList Ptg ----------------------------------
   *
   */
  PIOSTRING *Ptg;
  PIOBYTE flag_Ptg;  /* ==_PAR_TRUE if Ptg is present */
  PIOLONG n_Ptg; /* ==number_of_Ptg */

  /* ---- Optional InputList Hit ----------------------------------
   *
   */
  PIOSTRING *Hit;
  PIOBYTE flag_Hit;  /* ==_PAR_TRUE if Hit is present */
  PIOLONG n_Hit; /* ==number_of_Hit */

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

  /* ---- Optional OutputList Out_MAP --------------------------------------
   *
   */
  PIOSTRING *Out_MAP;
  PIOBYTE flag_Out_MAP;  /* ==_PAR_TRUE if Out_MAP is present */
  PIOLONG n_Out_MAP; /* ==number_of_Out_MAP */

  /* ---- Optional OutputList Out_VEC -------------------------------
   *
   */
  PIOSTRING *Out_VEC;
  PIOBYTE flag_Out_VEC;  /* ==_PAR_TRUE if Out_VEC is present */
  PIOLONG n_Out_VEC; /* ==number_of_Out_VEC */


  /* ---- Optional Parameter verbose -----------------------------------
   * Verbosity level. 0 : normal, 1:verbose
   */
  PIOINT verbose;
  PIOBYTE flag_verbose;  /* ==_PAR_TRUE if verbose is present */
  
  /* ---- Optional Parameter do_offset -----------------------------------
   * Verbosity level. 0 : normal, 1:do_offset
   */
  PIOINT do_offset;
  PIOBYTE flag_do_offset;  /* ==_PAR_TRUE if do_offset is present */

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
   * Parameter files for TOI simulations, one per bolo
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

  /* ---- Optional InputList addHPR_name -------------------------------
   * list of HPR object names of type PIOFLOAT to add to Signal input HPRs
   * must be a multiple of number of bolometers
   */
  PIOSTRING *addHPR_name;
  PIOBYTE   flag_addHPR_name;
  PIOLONG   n_addHPR_name;

  /* ---- Optional InputList addHPR_factor -------------------------------
   * list of multiplicative factor for each addHPR_name
   * must have the same number of elements as addHPR_name
   * if only one value is given, it will be used for all addHPR_name objects
   * default is 1.0 (simply add the addHPR_name to Signal)
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

  /* ---- Optional Parameter delta_psi ---------------------------------------------
   * angle in degrees to add to all bolometers psi pointing ({dbpath}/{pixname}_REP6_ptg_TUPLE_2)
   * to simulate an entire focal plane rotation.
   * omit or set to 0.0 for default/legacy behavior (no focal plane rotation)
   */
  PIOFLOAT *delta_psi;
  PIOBYTE flag_delta_psi;  /* = 1 if delta_psi is present */
  PIOLONG n_delta_psi;
  /* ---- Optional Parameter MAP_CNN_NETCDF ---------------------------------------------
  */
  PIOSTRING MAP_CNN;
  PIOBYTE flag_MAP_CNN;  /* ==_PAR_TRUE if MAP_CNN_NETCDF is present */
  /* ---- Optional Parameter INST_CNN_NETCDF ---------------------------------------------
  */
  PIOSTRING INST_CNN;
  PIOBYTE flag_INST_CNN;  /* ==_PAR_TRUE if INST_CNN_NETCDF is present */


} troll_parContent;

#endif
