#ifndef _STIM_PARAM_H_
#define _STIM_PARAM_H_


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>

#include "no_dmc_piolib_type_def.h"
#include "ring_index.h"
#include "brimo.h"

#include "RngStream.h" // http://arxiv.org/abs/1403.7645

#define ADCBITS  (16l)
#define INLSIZE  (1l<<ADCBITS)
#define MAXRSTREAMS (26100l)

/*******************************************************************************
parcontent: structure containing the read parameter file
WARNING!! this must be the first structure of the .h file, or fillParLoader.py won't work
*******************************************************************************/

typedef struct {
  PIOSTRING bolometer;
  PIOSTRING brimo_filename;
  PIOINT    stay_in_memory;
  PIOINT    iterations;
  PIOINT    BeginRing;
  PIOINT    EndRing;
  PIOINT    rings_per_read;

  // signal_in: name of signal to read from disk, set to any constant float value to work with a constant input signal
  PIOSTRING signal_in;

  // signal_out: name of signal to write to disk, set to '0' to not write output signal to disk
  PIOSTRING signal_out;

  // hpr_out: name of hpr to write to disk, set to '0' to not write output HPR to disk
  PIOSTRING hpr_out;

  // do_SHPR: set to 1 to produce splined HPR (APR20 and SRoll3)
  // or set to 0 to produce "legacy" HPR (JAN18 and SRoll1 and SRoll2.x)
  PIOINT do_SHPR;

  // pixel_number: TOI of pixel index in HPR for projection
  PIOSTRING TOI_pixel_index;

  // must contain the Phase_dx11 object name if despike_flag=1
  PIOSTRING phase;

  // random_seed is used for reproducibility of various stim modules using random generated numbers
  PIOINT    random_seed;

  // set to 0 to not add photonic noise, or 1 to add noise with level from BRIMO, or a float to override BRIMO white noise value (KCMB)
  PIOFLOAT  do_photonic_noise;

  PIOINT    do_LFER_convolve;

  // set nharm_oof to 0 to not add one over f noise, or N to add oof noise and give N amplitudes for the first N harmonics
  PIOINT    nharm_oof;
  PIOSTRING amp_oof_param;

  // set to 0 to not add electronic noise, or 1 to add noise with level from BRIMO, or a float to override BRIMO white noise value (KCMB)
  PIOFLOAT  do_electronic_noise;

  PIOFLOAT  do_despike_flag;

  // set to 1 to convert signal to ADU, add baseline and raw constant, modulate and add ADC nonlinearity effect
  PIOINT    do_sim_adu;
  PIOSTRING sim_inl_name;
  PIOSTRING rawgain_name;
  PIOSTRING rawcst_name;
  PIOSTRING fourk_name;
  PIOSTRING add_baseline;
  PIOFLOAT  electronic_noise_adu;

  PIOINT    do_compress_decompress;

  PIOINT    do_correct_adc;
  PIOSTRING corr_inl_name;

//  // if not 0, contains an object name which content will be subtracted from corrected ADU
//  PIOSTRING do_subtract_adu;

  PIOINT    do_adu_to_volts;
  
  PIOINT    do_bl_demod;

  PIOINT    do_gaindecorr;
  PIOINT    use_bolometer_nonlinearity;
  PIOSTRING thermal_baseline;

  PIOINT    do_LFER_deconvolve;
  PIOSTRING deconv_LPF_FILTER;  // set to "NONE" to not use the lowpass filter at deconvolution for freqs < 545, set to '0' to use default value
  PIOSTRING deconv_R_FILTER;    // set to "NONE" to not use the regularisation filter at deconv for freqs >= 545, set to '0' to use default value

  // if set to 1, will convert signal from Watts to KCMB (or MJY/sr at 545/857) using BRIMO calibration factor
  // if set to -1, will convert signal fro KCMB (or MJY/sr at 545/857) to Watts
  // SRoll input HPR must be in Watts.
  // stim "standard" output signal unit is Watts at 100-353 and MJy/sr at 545/857 (545/857 don't go through ADU/ADC)
  // so to produce Watts HPR for SRoll with standard simulations, do_calibrate must be set to 0 for 100-353 and -1 for 545/857
  // if set to neither 0, 1 or -1, will multiply the signal by this value before saving TOI and converting to HPR
  PIOFLOAT  do_calibrate;

  // if set, the TOI in add_final_toi will be multiplied by add_final_toi_factor and added to
  // simulated TOI after calibration and before projection to HPR
  PIOSTRING add_final_toi;
  PIOFLOAT  add_final_toi_factor;

  // optional list of HPRs to add to the final stim HPR returned to sroll
  PIOSTRING *add_hpr_name;
  PIOBYTE flag_add_hpr_name;  // == 1 if add_hpr_name is present, else 0
  PIOLONG n_add_hpr_name;     // == number_of_add_hpr_name

  PIOFLOAT *add_hpr_factor;
  PIOBYTE flag_add_hpr_factor;  // == 1 if add_hpr_factor is present, else 0
  PIOLONG n_add_hpr_factor;     // == number_of_add_hpr_factor

  PIOINT  RD12calib;       // set to 1 to use RD12 calib factors (>=2018), or 0 or omitted for FFP10 reproduction
  PIOBYTE flag_RD12calib;  // == 1 if RD12calib is present, else 0

  PIOSTRING save_sim_inl;       // set to a file name prefix to save the simulation INL to this binary file, will be postfixed with "_{pixname}_{iternum:%03d}.float64"
  PIOBYTE   flag_save_sim_inl;  // == 1 if save_sim_inl is present, else 0

} stim_parContent;


/*******************************************************************************
stim_Data: structure to hold in memory all data used by stim
to iterate without touching disks
*******************************************************************************/

typedef struct {
  PIOFLOAT      *signal_in;        // 4 Bytes per sample
  PIOINT        *pixel_index;      // 4 Bytes per sample
  PIODOUBLE     *phase;            // 8 Bytes per sample
  PIODOUBLE     *sim_inl;          // 512KB
  int           raw_begin_ring;    // 4 Bytes
  PIODOUBLE     *raw_gain;         // 640 Bytes
  PIODOUBLE     *raw_cst;          // 640 Bytes per ring
  int           n_harmonics;       // 4 Bytes
  float         *fourk_harmonics;  // 48 Bytes
  float complex *fourk_amplitudes; // 96 Bytes per ring
  PIOFLOAT      *dsn_baseline;     // 4 Bytes per sample
  PIODOUBLE     *corr_inl;         // 512KB
  PIOFLOAT      *thermal_baseline; // 4 Bytes per sample
  PIOFLOAT      *add_toi;          // 4 Bytes per sample
  PIOFLOAT      **add_hpr;         // 108KB per ring per hpr
} stim_Data;


/*******************************************************************************
stimParameters: structure containing all necessary data for adding and removing
effects modules called by stim.c.
*******************************************************************************/

typedef struct {
  stim_parContent Param;                 // parameter file content
  BRIMO      brimo;                      // binary reduced IMO values
  stim_Data  data;                       // input data cache
  int        mpi_size, mpi_rank;
  PIOSTRING  msg_prefix;                 // prefix that can be used in messages to know who's speaking
  int        trace_level;                // 0: no trace, 1: tracing function calls, 2: more details
  RngStream  rstream[MAXRSTREAMS];       // one independant RNG stream per ring
  long       signal_first_sample_number; // buffer first sample number
  long       total_signal_length;        // length of buffer
  int        BeginRing;                  // first ring number to process
  int        EndRing;                    // last first ring number to process
  long       beginring_first_sample;     // BeginRing first sample offset in buffer
  long       endring_last_sample;        // EndRing last sample offset in buffer
  long long  local_seed;                 // random seed initialised by stim for each effect
  double     calib;                      // calibration factor (KCMB2WATT), either RD12 (>=2018) or DX11 (FFP10) values
  int        conv_CONVOLVE;   // convolve.c convolution direction
  PIOSTRING  conv_R_FILTER;   // convolve.c regularisation filter
  PIOSTRING  conv_LPF_FILTER; // convolve.c low-pass filter
  long       conv_start;      // sample number where convolution starts
  long       conv_end;        // sample number where convolution ends
  long       deconv_start;    // sample number where deconvolution starts
  long       deconv_end;      // sample number where deconvolution ends
  long       first_full_ring; // first full ring to process in despike_flag, sim_adu and correct_adc
  long       last_full_ring;  // last full ring to process in despike_flag, sim_adu and correct_adc
  long       smooth_start;    // sample number where baseline smoothing starts
  long       smooth_end;      // sample number where baseline smoothing ends
} stimParameters;

#endif
