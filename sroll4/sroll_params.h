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
 
#ifndef sroll_PARAM_H_
#define sroll_PARAM_H_


#include "sroll_param.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


/*
 * parcontent : structure containing the read parameter file
 */

typedef struct {
	PIOINT number_of_Out_xi2_corr;
	PIOSTRING *Out_xi2_corr;
	PIOINT number_of_Out_xi2;
	PIOSTRING *Out_xi2;
	PIOINT number_of_Out_Offset_corr;
	PIOSTRING *Out_Offset_corr;
	PIOINT number_of_Out_Offset;
	PIOSTRING *Out_Offset;
	PIOINT n_MAP;
	PIOINT number_of_MAP;
	PIOSTRING *MAP;
	PIOSTRING *bolo_map;
	PIOINT number_of_rgcnn;
	PIOSTRING *rgcnn;
	PIOINT number_of_phase;
	PIOSTRING *phase;
	PIOINT number_of_fsl;
	PIOSTRING *fsl;
	PIOINT number_of_Theo_noPS;
	PIOSTRING *Theo_noPS;
	PIOINT number_of_Ptg_noPS;
	PIOSTRING *Ptg_noPS;
	PIOINT number_of_Hit_noPS;
	PIOSTRING *Hit_noPS;
	PIOINT number_of_DipOrb_noPS;
	PIOSTRING *DipOrb_noPS;
	PIOINT number_of_Badring;
	PIOSTRING *Badring;
	PIOINT number_of_ADU;
	PIOSTRING *ADU;
	PIOINT number_of_Signal_noPS;
	PIOSTRING *Signal_noPS;
	PIOINT number_of_in_polar_fit_U;
	PIOSTRING *in_polar_fit_U;
	PIOINT number_of_in_polar_fit_Q;
	PIOSTRING *in_polar_fit_Q;
	PIOINT number_of_in_template_map_U;
	PIOSTRING *in_template_map_U;
	PIOINT number_of_in_template_map_Q;
	PIOSTRING *in_template_map_Q;
	PIOINT number_of_in_template_map_I;
	PIOSTRING *in_template_map_I;
	PIOSTRING TEMPLATEMAP;
	PIOSTRING Mask;
	PIOINT number_of_MAPRINGS;
	PIOINT *MAPRINGS;
	PIOINT number_of_bolomask;
	PIOINT *bolomask;
	PIOINT n_OUT_NOPOL;
	PIOINT number_of_NOPOL;
	PIOINT *OUT_NOPOL;
	PIOINT number_of_FSLCOEF;
	PIOFLOAT *FSLCOEF;
	PIOINT number_of_Monop;
	PIOFLOAT *Monop;
	PIOINT number_of_NEP;
	PIOFLOAT *NEP;
	PIOINT number_of_CrossPol;
	PIOFLOAT *CrossPol;
	PIOINT number_of_Calibration;
	PIOFLOAT *Calibration;
	PIOSTRING MAP_CNN;
	PIOSTRING INST_CNN;
	PIOINT N_IN_ITT;
	PIOFLOAT CNN_RESIDU;
	PIOINT CNN_YSIZE;
	PIOINT CNN_XSIZE;
	PIOSTRING CALLCNN;
	PIOINT number_of_DOCNN;
	PIOINT *DOCNN;
	PIOINT *NADUSTEP;
	PIOINT number_of_NADUSTEP;
	PIOINT *NADU;
	PIOINT number_of_NADU;
	PIOINT GAINSTEP;
	PIOINT number_of_ADDCO13;
	PIOINT SEED1;
	PIOINT number_of_SEED;
	PIOINT BUILDTF;
	PIOINT TEMPLATE_NSIDE;
	PIOINT verbose;
	PIOFLOAT delta_psi;
	PIOINT saveCO;
	PIOINT saveCOV;
	PIOINT FITPOLEFF;
	PIOINT FITANGLE;
	PIOINT NITT;
	PIOFLOAT seuilcond;
	PIOFLOAT XI2STOP;
	PIOINT TESTPOL;
	PIOINT SAVEINTMAP;
	PIOINT KCMBIN;
	PIOINT ADDDIP;
	PIOINT REMHDIP;
	PIOINT REMDIP;
	PIOINT Nside;
	PIOINT D_NOPOL;
	PIOINT DOMAXVRAIE;
	PIOINT DOGAINDIP;
	PIOINT DODIPCAL;
	PIOINT CUTRG;
	PIOINT CALCODUST;
	PIOINT AVG13CO;
	PIOINT AVG12CO;
	PIOINT AVFCO;
	PIOINT AVGDUST;
	PIOINT AVFDUST;
	PIOINT AVGDUST100;
	PIOINT AVGSYNCHRO;
	PIOINT AVFSYNC;
	PIOINT AVGFREEFREE;
	PIOINT AVFFF;
	PIOINT RSTEP;
	PIOINT EndRing;
	PIOINT BeginRing;
	PIOSTRING plop;
	PIOSTRING *bolo2;
	PIOSTRING *bolo;
	PIOINT nbolo;
}sroll_parContent

#endif
