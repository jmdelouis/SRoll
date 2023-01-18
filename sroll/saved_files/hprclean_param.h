/*
 *
 * Template for hprclean parameter file reader include file
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

#ifndef _HPRCLEAN_PARAM_H_
#define _HPRCLEAN_PARAM_H_


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

  /* ---- Input BEGINRINGINDEX ---------------------------------------
   * Input BEGINRINGINDEX index of the input timeline
   */

  PIOSTRING BEGINRINGINDEX;
  /* ---- InputList HIDX ---------------------------------------
   * Input HIDX index to recompute HPR
   * one per bolometer, and one bolometer per MPI rank
   */
  PIOSTRING *HIDX;
  PIOINT n_HIDX;

  /* ---- ParameterList Calib ---------------------------------------
   * ParameterList
   * Same calibration factors from KCMB to Watts as used in stim and SRoll,
   * used to convert zodi from KCMB to Watts to subtract it from Signal in Watts
   * one per bolometer, and one bolometer per MPI rank
   */
  PIODOUBLE *calib;
  PIOINT n_calib;

  /* ---- InputList Zodi ---------------------------------------
   * InputList Zodi model to remove
   * one per bolometer, and one bolometer per MPI rank
   */
  PIOSTRING *zodi;
  PIOBYTE flag_zodi;  /* == 1 if zodi is present */
  PIOINT n_zodi;

  /* ---- InputList Signal ---------------------------------------
   * InputList Signal 
   * one per bolometer, and one bolometer per MPI rank
   */
  PIOSTRING *Signal;
  PIOINT n_Signal;

  /* ---- OutputList hpr ---------------------------------------
   * OutputList hpr (should be equivalent to REP6)
   * one per bolometer, and one bolometer per MPI rank
   */
  PIOSTRING *hpr;
  PIOBYTE flag_hpr;  /* == 1 if hpr is present */
  PIOINT n_hpr;

  /* ---- OutputList hpr_corr ---------------------------------------
   * OutputList hpr clean from large scale noise (based on spline destriping)
   * one per bolometer, and one bolometer per MPI rank
   */
  PIOSTRING *hpr_corr;
  PIOINT n_hpr_corr;


} hprclean_parContent;

#endif
