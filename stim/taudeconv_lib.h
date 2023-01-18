
#ifndef _TAUDECONVLIB_H
#define _TAUDECONVLIB_H


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "no_dmc_piolib_type_def.h"


// defining tau deconvolver parameters within a single structure
struct par_taudeconv{

  char method[256];
  char r_method[256];
  char n_method[256];
  char lpf_method[256];

// ------- Convolution flag
  PIOBYTE   convol; // set to 1 for convolution, 0 for deconvolution

// --- General parameters
  PIODOUBLE fsampling;

// ---- For Time constant method
  PIODOUBLE tau1, tau2, tau3, tau4, tau5, tau6, tau7, tau8;
  PIODOUBLE A1, A2, A3, A4, A5, A6, A7, A8;
  PIODOUBLE res1, res2, res3, ele1, ele2, ele3;
  PIODOUBLE tau_stray, sphase,tauhp;
  PIODOUBLE tfmodclip;

// ---- For Input Filter
  PIOLONG   nbDataVect;
  PIODOUBLE *freq, *amp, *phase;  

// --- time delay parameters in units of samples
  PIODOUBLE sphase_offset, global_offset, total_offset; 

// ---- Paarameters for Regularization filtering -----
  PIODOUBLE fwidth_rfilter;
  
// ----- Parameters for notch filter ------
  PIOLONG    nfilter_nfreq;
  PIODOUBLE *nfilter_fval, *nfilter_fwidth;

// ---- Paarameters for LPF filter  -----
// ----------- EXP * COS^2 filter -----------
  PIODOUBLE fc_lpf, fgauss_lpf, ffact_lpf;
  
};

// Function interfaces

void array_complex_mult(PIOLONG Ndata, PIODOUBLE *toi_1, PIODOUBLE *toi_2);
void construct_timeresponse_inv(struct par_taudeconv pardeconv, PIOLONG Ndata, PIODOUBLE *tresponse);
void construct_regularization_filter(struct par_taudeconv pardeconv, PIOLONG Ndata, PIODOUBLE *tresponse);
void construct_notch_filter(struct par_taudeconv pardeconv, PIOLONG Ndata, PIODOUBLE *tresponse);
void construct_lpf_filter(struct par_taudeconv pardeconv, PIOLONG Ndata, PIODOUBLE *tresponse);
void tf_resonance_JH(struct par_taudeconv pardeconv, PIOLONG Ndata, PIODOUBLE *fr, PIODOUBLE *toi);
void tf_electronic_JH(struct par_taudeconv pardeconv, PIOLONG Ndata, PIODOUBLE *fr, PIODOUBLE *toi);
void lowpass(PIODOUBLE w, PIODOUBLE tau, PIODOUBLE *zr, PIODOUBLE *zi);
void set_complex_exp(PIODOUBLE mod, PIODOUBLE arg,PIODOUBLE *zr,PIODOUBLE *zi);
void complex_multiplication(PIODOUBLE *r1, PIODOUBLE *i1, PIODOUBLE r2, PIODOUBLE i2);
void complex_division(PIODOUBLE *r1, PIODOUBLE *i1, PIODOUBLE r2, PIODOUBLE i2);
void tf_modclipping(struct par_taudeconv pardeconv, PIOLONG Ndata, PIODOUBLE *toi); 
void tf_write(PIOSTRING filename, PIOLONG Ndata, PIODOUBLE *fr, PIODOUBLE *toi);
void tf_8time_constants_fourier(struct par_taudeconv pardeconv, PIOLONG Ndata, PIODOUBLE *fr, PIODOUBLE *toi);
// Added by B. Crill 8July2010
void tf_jh8(struct par_taudeconv pardeconv,PIOLONG Ndata, PIODOUBLE *fr, PIODOUBLE *toi);
void tf_elect_anal_jh8(PIODOUBLE f, PIODOUBLE tau0, PIODOUBLE fmod, PIODOUBLE sphase, PIODOUBLE *tfr, PIODOUBLE *tfi);

/* removed fy mottet@iap.fr for stim

void tau_deconv(stim_parContent *Param);
void timeresp_deconv(struct par_taudeconv pardeconv, PIOLONG Ndata, PIODOUBLE *toideconv);
void read_imofile( stim_parContent *Param, struct par_taudeconv *tauparams);
void interp_invert_filter(PIOLONG Ndata, PIODOUBLE *fr,PIOLONG nbDataVect, PIODOUBLE *freq, PIODOUBLE *amp,PIODOUBLE *phase, PIODOUBLE *tresponse);
void tf_3time_constants(struct par_taudeconv pardeconv, PIOLONG Ndata,PIODOUBLE *time, PIODOUBLE *tresponse);
void tf_3time_constants_fourier(struct par_taudeconv pardeconv, PIOLONG Ndata,PIODOUBLE *fr, PIODOUBLE *toi);
void tf_4time_constants_fourier(struct par_taudeconv pardeconv, PIOLONG Ndata,PIODOUBLE *fr, PIODOUBLE *toi);
void tf_5time_constants_fourier(struct par_taudeconv pardeconv, PIOLONG Ndata,PIODOUBLE *fr, PIODOUBLE *toi);
void tf_6time_constants_fourier(struct par_taudeconv pardeconv, PIOLONG Ndata,PIODOUBLE *fr, PIODOUBLE *toi);

void tf_roma_analytic(struct par_taudeconv pardeconv, PIOLONG Ndata, PIODOUBLE *fr, PIODOUBLE *toi);
void tf_elect_anal_roma(PIODOUBLE f, PIODOUBLE tau0, PIODOUBLE tauhp, PIODOUBLE fmod, 
                        PIODOUBLE sphase,PIODOUBLE *tfr, PIODOUBLE *tfi);
void hfiele_grenoble(PIODOUBLE w, PIODOUBLE tau_stray, PIODOUBLE tauhp, PIODOUBLE *rval, PIODOUBLE *ival);

*/

#endif
