
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <sys/time.h>
#include <unistd.h>   
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "no_dmc_metadata.h"
#include "no_dmc_data_access.h"
#include "no_dmc_piolib_type_def.h"
#include "no_dmc_util.h"
#include "no_dmc_debug.h"
#include "no_dmc_version.h"

#include "stim_param.h"
#include "stim_parLoader.h"
#include "binary_file_toi.h"
#include "ring_index.h"
#include "stim_tools.h"

#include "add_oof_and_white_noise.h"
#include "add_oof_noise.h"
#include "real_filter.h"
#include "convolve.h"
#include "sim_adu.h"
#include "compress_decompress.h"
#include "correct_adc.h"
#include "adu_to_volts.h"
#include "despike_flag.h"
#include "bl_demod.h"
#include "gain_decorr.h"

#include "../sroll/HPR_Cleaner.h"

#define MAXSAMPLEPERRING (785000l)


/*******************************************************************************
init_stim()
function reading stim parameter file and brimo and common input data
and broadcast them to all ranks
return 0 when ok
*******************************************************************************/

int init_stim( stimParameters *stimPar, PIOSTRING param_filename) {

  stim_parContent *Param = &stimPar->Param;
  stim_Data       *data  = &stimPar->data;
  BRIMO           *brimo = &stimPar->brimo;

  int i, mpi_rank, mpi_size;

  MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size( MPI_COMM_WORLD, &mpi_size);

  // read parameter file and brimo only once with rank 0
  if (mpi_rank == 0) {

    // read parameters
    if (stim_readParam( &stimPar->Param, param_filename) != 0) {
      printf( "%s: unable to parse the parameter file %s.\n", __FILE__, param_filename);
      exit( -1);
    }
    assert( Param->n_add_hpr_name == Param->n_add_hpr_factor); // just to check...
    readBRIMO( Param->brimo_filename, getBC( Param->bolometer), &stimPar->brimo);

    // get a trace of brimo values used
    printBRIMO( &stimPar->brimo);

    // select calibration
    if ((Param->flag_RD12calib == 1) && (Param->RD12calib == 1)) {
      // simulations from 2018 and on use RD12 calibration factors
      stimPar->calib = get_RD12calib( brimo->pixname);
      printf( "%s: using RD12 calib=%g.\n", __FILE__, stimPar->calib);
    }
    else {
      // FFP10 and "fiducial/dustapr17" 2016-2017 simulations use DX11 calibration factors
      stimPar->calib = brimo->KCMB2WATT;
      printf( "%s: using DX11 calib=%g.\n", __FILE__, stimPar->calib);
    }

    // initialise data pointers
    data->signal_in = NULL;
    data->pixel_index = NULL;
    data->phase = NULL;
    data->sim_inl = NULL;
    data->raw_begin_ring = 0; // starting ring of raw data in memory, useful while sim_adu and correct_adc use same inputs but not on same ring range
    data->raw_gain = NULL;
    data->raw_cst = NULL;
    data->n_harmonics = 0;
    data->fourk_harmonics = NULL;
    data->fourk_amplitudes = NULL;
    data->dsn_baseline = NULL;
    data->corr_inl = NULL;
    data->thermal_baseline = NULL;
    data->add_toi = NULL;
    data->add_hpr = NULL;
  }

  // broadcast parameter file and brimo to all ranks
  assert( MPI_Bcast( stimPar, sizeof( *stimPar), MPI_BYTE, 0, MPI_COMM_WORLD) == MPI_SUCCESS);

  // set local mpi_rank and size
  stimPar->mpi_size = mpi_size;
  stimPar->mpi_rank = mpi_rank;
  int ndec = (int) ceil( log10( mpi_size));
  if (ndec == 0) ndec=1;
  sprintf( stimPar->msg_prefix, "rank %.*d/%d:", ndec, mpi_rank, mpi_size);
  if (mpi_rank == 0) {
    stimPar->trace_level = 1;
  } else {
    stimPar->trace_level = 0;
  }

  // broadcast add_hpr* lists
  if (Param->flag_add_hpr_name == _PAR_TRUE) {
    assert( (Param->flag_add_hpr_factor == _PAR_TRUE)
         && (Param->n_add_hpr_factor == Param->n_add_hpr_name));
    if (mpi_rank != 0) {
      Param->add_hpr_name = malloc(Param->n_add_hpr_name * sizeof( PIOSTRING));
      assert( Param->add_hpr_name != NULL);
      Param->add_hpr_factor = malloc(Param->n_add_hpr_factor * sizeof( PIOFLOAT));
      assert( Param->add_hpr_factor != NULL);
    }
    assert( MPI_Bcast( Param->add_hpr_name,   Param->n_add_hpr_name * sizeof( PIOSTRING), MPI_BYTE, 0, MPI_COMM_WORLD) == MPI_SUCCESS);
    assert( MPI_Bcast( Param->add_hpr_factor, Param->n_add_hpr_name * sizeof( PIOFLOAT),  MPI_BYTE, 0, MPI_COMM_WORLD) == MPI_SUCCESS);
  }

  // read and broadcast sim_inl
  if (Param->do_sim_adu == 1) {
    if (strcmp( Param->sim_inl_name, "0") != 0) {
      data->sim_inl = malloc( INLSIZE * sizeof( PIODOUBLE));
      if (mpi_rank == 0) {
        assert( noDMC_readObject_PIODOUBLE( Param->sim_inl_name, 0, INLSIZE, data->sim_inl) >= 0);
      }
      assert( MPI_Bcast( data->sim_inl, INLSIZE * sizeof( PIODOUBLE), MPI_BYTE, 0, MPI_COMM_WORLD) == MPI_SUCCESS);
    }
  }

  // read, process and broadcast raw_gain
  if (Param->do_sim_adu == 1) {
    PIODOUBLE tmp_gain[80];
    data->raw_gain = malloc( 80 * sizeof( PIODOUBLE));
    if (mpi_rank == 0) {
      assert( noDMC_readObject_PIODOUBLE( Param->rawgain_name, 0, 80, tmp_gain) >= 0);
      double gain_norm = 0.0;
      for (i = 0; i < 80; i++) {
        data->raw_gain[i] = tmp_gain[(i + brimo->sphase) % 80];
      }
      for (i = 0; i < 40; i++) gain_norm += data->raw_gain[i];
      gain_norm = fabs( gain_norm);
      for (i = 0; i < 80; i++) {
        data->raw_gain[i] /= gain_norm;
      }  
    }
    assert( MPI_Bcast( data->raw_gain, 80 * sizeof( PIODOUBLE), MPI_BYTE, 0, MPI_COMM_WORLD) == MPI_SUCCESS);
  }

  // read and broadcast corr_inl
  if ((Param->do_sim_adu == 1) || (Param->do_correct_adc == 1)) {
    data->corr_inl = malloc( INLSIZE * sizeof( PIODOUBLE));
    if (strcmp( Param->corr_inl_name, "0") != 0) {
      if (mpi_rank == 0) {
        assert( noDMC_readObject_PIODOUBLE( Param->corr_inl_name, 0, INLSIZE, data->corr_inl) >= 0);
      }
      assert( MPI_Bcast( data->corr_inl, INLSIZE * sizeof( PIODOUBLE), MPI_BYTE, 0, MPI_COMM_WORLD) == MPI_SUCCESS);
    } else {
      for (i = 0; i < INLSIZE; i++) {
        data->corr_inl[i] = i;
      }
    }
  }

  return( 0);
}


/*******************************************************************************
free_stim()
release all memory allocated by stim
*******************************************************************************/

int free_stim( stimParameters *stimPar) {

  stim_parContent *Param = &stimPar->Param;
  stim_Data       *data  = &stimPar->data;
  int i;

  if (Param->stay_in_memory) {
    STIM_TRACE( 1, "freeing stim data");
    if (data->signal_in != NULL) {
      free( data->signal_in);
      data->signal_in = NULL;
    }
    if (data->pixel_index != NULL) {
      free( data->pixel_index);
      data->pixel_index = NULL;
    }
    if (data->phase != NULL) {
      free( data->phase);
      data->phase = NULL;
    }
    if (data->sim_inl != NULL) {
      free( data->sim_inl);
      data->sim_inl = NULL;
    }
    if (data->raw_gain != NULL) {
      free( data->raw_gain);
      data->raw_gain = NULL;
    }
    if (data->raw_cst != NULL) {
      free( data->raw_cst);
      data->raw_cst = NULL;
    }
    data->raw_begin_ring = 0;
    if (data->fourk_harmonics != NULL) {
      free( data->fourk_harmonics);
      data->fourk_harmonics = NULL;
    }
    if (data->fourk_amplitudes != NULL) {
      free( data->fourk_amplitudes);
      data->fourk_amplitudes = NULL;
    }
    if (data->dsn_baseline != NULL) {
      free( data->dsn_baseline);
      data->dsn_baseline = NULL;
    }
    if (data->corr_inl != NULL) {
      free( data->corr_inl);
      data->corr_inl = NULL;
    }
    if (data->thermal_baseline != NULL) {
      free( data->thermal_baseline);
      data->thermal_baseline = NULL;
    }
    if (data->add_toi != NULL) {
      free( data->add_toi);
      data->add_toi = NULL;
    }
    if (data->add_hpr != NULL) {
      for (i = 0; i < Param->n_add_hpr_name; i++) {
        free( data->add_hpr[i]);
      }
      free( data->add_hpr);
      data->add_hpr = NULL;
    }
  }

  if (Param->flag_add_hpr_name == _PAR_TRUE) {
    free( Param->add_hpr_name);
    Param->add_hpr_name = NULL;
    free( Param->add_hpr_factor);
    Param->add_hpr_factor = NULL;
    free( data->add_hpr);
    data->add_hpr = NULL;
  }

  return( 0);
}


/*******************************************************************************
stim()
function calling all functions adding and removing effects
return 0 when ok
*******************************************************************************/

int stim( stimParameters *stimPar, int begin_ring, int rings_to_process, PIOFLOAT *hpr_out) {

  stim_parContent *Param = &stimPar->Param;
  stim_Data       *data  = &stimPar->data;
  BRIMO           *brimo = &stimPar->brimo;
  long   i;
  struct timeval  t0;
  PIOSTRING outname;

  if (stimPar->trace_level >= 1) {
    gettimeofday( &t0, NULL);
    fprintf( stderr, "\n");
    STIM_TRACE( 1, "%s, seeed=%d, processing rings (%d-%d) (%d rings, %de6 samples)", Param->bolometer, Param->random_seed, begin_ring, begin_ring + rings_to_process-1, rings_to_process, (int)((ENDRINGINDEX(begin_ring + rings_to_process-1)-BEGINRINGINDEX(begin_ring)+1)/1e6));
  }

  assert( hpr_out != NULL);

  // if set to auto-value (-1) and not set by sroll, don't keep inputs in memory
  if (Param->stay_in_memory == -1) {
    Param->stay_in_memory = 0;
  }

  // test if signal_in is an object name or a constant float value
  char *endptr;
  float constant_signal = strtod( Param->signal_in, &endptr);
  if (endptr != Param->signal_in) {
    if (stimPar->mpi_rank == 0) {
      STIM_TRACE( 0, "detected float value (%g) in Param->signal_in, using a constant as input", constant_signal);
    }
    Param->signal_in[0] = 0;
  }


//==============================================================================
// initialise RngStream with iteration number and bolometer BC

  // +1 because seeds can't be 0
  unsigned long germe[6] = { Param->random_seed+1, brimo->bc+1,
                             Param->random_seed+1, brimo->bc+1,
                             Param->random_seed+1, brimo->bc+1 };
  RngStream_SetPackageSeed( germe);

//==============================================================================
// loop on rings by Param->rings_per_read increments

  int rings_per_read;
  if (Param->rings_per_read == -1) {
    rings_per_read = rings_to_process;
  } else {
    rings_per_read = Param->rings_per_read;
  }

  for (int loop_begin_ring = begin_ring; loop_begin_ring < begin_ring + rings_to_process; loop_begin_ring += rings_per_read) {
    // first and last sample to process, will be written in output
    int loop_end_ring = loop_begin_ring + rings_per_read - 1;
    if (loop_end_ring > begin_ring + rings_to_process - 1) {
      loop_end_ring = begin_ring + rings_to_process - 1;
    }
    long loop_first_sample = BEGINRINGINDEX( loop_begin_ring);
    long loop_last_sample  = ENDRINGINDEX(   loop_end_ring);
    long loop_sample_count = loop_last_sample - loop_first_sample + 1;


//==============================================================================
// compute margins to correctly process ring range

    long read_first_sample = loop_first_sample;
    long read_last_sample  = loop_last_sample;

    if (Param->do_bl_demod == 1) {
      // baseline demodulation needs a LSMOOTH/2 margin around sample interval to process
      stimPar->smooth_start = read_first_sample;
      stimPar->smooth_end   = read_last_sample;
      read_first_sample -= LSMOOTH / 2;
      read_last_sample  += LSMOOTH / 2;
    }

    if (Param->do_sim_adu == 1 || Param->do_correct_adc == 1) {
      // find the interval of full rings to include read_first/last_sample
      int br = loop_begin_ring;
      while (BEGINRINGINDEX( br) > read_first_sample) br--;
      read_first_sample = BEGINRINGINDEX( br);
      int er = loop_end_ring;
      while (ENDRINGINDEX( er) < read_last_sample) er++;
      read_last_sample = ENDRINGINDEX( er);
      stimPar->first_full_ring = br;
      stimPar->last_full_ring  = er;
    }

    if (Param->do_LFER_deconvolve == 1) {
      // add margins before and after read_first/last_sample to correctly process the loop interval
      long first_conv_start = read_first_sample / LDECONV * LDECONV;
      int convolution_count = (read_last_sample - first_conv_start + 1) / LDECONV + 1;
      long margin_first_sample = first_conv_start - HALFPERIOD;
      long margin_last_sample  = first_conv_start + convolution_count * LDECONV + HALFPERIOD - 1;
      // check that margins are correctly computed
      assert( (margin_last_sample - read_last_sample) >= HALFPERIOD);
      assert( (margin_last_sample - read_last_sample) < LDECONV + HALFPERIOD);
      assert( (margin_last_sample - margin_first_sample + 1) % LDECONV == 0);
      // update read interval
      read_first_sample = margin_first_sample;
      read_last_sample  = margin_last_sample;
      stimPar->deconv_start = first_conv_start;
      stimPar->deconv_end   = first_conv_start + convolution_count * LDECONV - 1;
    }

    if (Param->do_despike_flag != 0) {
      // find the interval of full rings to include read_first/last_sample
      int br = loop_begin_ring;
      while (BEGINRINGINDEX( br) > read_first_sample) br--;
      read_first_sample = BEGINRINGINDEX( br);
      int er = loop_end_ring;
      while (ENDRINGINDEX( er) < read_last_sample) er++;
      read_last_sample = ENDRINGINDEX( er);
      stimPar->first_full_ring = br;
      stimPar->last_full_ring  = er;
    }

    if (Param->do_LFER_convolve == 1) {
      // add margins before and after read_first/last_sample to correctly process the loop interval
      long first_conv_start = read_first_sample / LDECONV * LDECONV;
      int convolution_count = (read_last_sample - first_conv_start + 1) / LDECONV + 1;
      long margin_first_sample = first_conv_start - HALFPERIOD;
      long margin_last_sample  = first_conv_start + convolution_count * LDECONV + HALFPERIOD - 1;
      // check that margins are correctly computed
      assert( (margin_last_sample - read_last_sample) >= HALFPERIOD);
      assert( (margin_last_sample - read_last_sample) < LDECONV + HALFPERIOD);
      assert( (margin_last_sample - margin_first_sample + 1) % LDECONV == 0);
      // update read interval
      read_first_sample = margin_first_sample;
      read_last_sample  = margin_last_sample;
      stimPar->conv_start = first_conv_start;
      stimPar->conv_end   = first_conv_start + convolution_count * LDECONV - 1;
}

#if 0
    // print data range and margins per rank information
    for (int irank=0; irank<stimPar->mpi_size; irank++) {
      if (stimPar->mpi_rank == irank) {
        printf( "%s ring range: %d-%d, %d rings, %.2fe6 samples, %.2fe6 margin samples\n",
                stimPar->msg_prefix,
                loop_begin_ring, loop_end_ring,
                loop_end_ring - loop_begin_ring + 1,
                loop_sample_count / 1e6,
                (read_last_sample - read_first_sample + 1 - loop_sample_count) / 1e6);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
#endif



//==============================================================================
// read input data

    // read input signal if not already in memory, else reuse read data
    PIOLONG read_sample_count = read_last_sample - read_first_sample + 1;
    PIOFLOAT *signal = NULL;
    if ((!Param->stay_in_memory) || (data->signal_in == NULL)) {
      STIM_TRACE( 1, "reading %ld samples before and %ld samples after processing interval", loop_first_sample - read_first_sample, read_last_sample - loop_last_sample);
      signal = malloc( read_sample_count * sizeof( PIOFLOAT));
      assert( (signal != NULL) && "Not enough memory to allocate signal array");
      if (Param->signal_in[0] != 0) {
        assert( noDMC_readObject_PIOINTorPIOFLOAT( Param->signal_in, read_first_sample, read_sample_count, signal) == 0);
      } else {
        for (i = 0; i < read_sample_count; i++) {
          signal[i] = constant_signal;
        }
      }
    }
    if (Param->stay_in_memory) {
      if (data->signal_in == NULL) {
        // first iteration, store signal in data struct
        data->signal_in = signal;
        signal = NULL;
      } else {
        STIM_TRACE( 1, "reusing input data kept in memory");
      }
      // stay in memory, (re)use data.signal_in
      signal = malloc( read_sample_count * sizeof( PIOFLOAT));
      assert( (signal != NULL) && "Not enough memory to allocate signal array");
      memcpy( signal, data->signal_in, read_sample_count * sizeof( PIOFLOAT));
    }

    // read TOI_pixel_index from loop_begin_ring to loop_end_ring
    if (data->pixel_index == NULL) {
      data->pixel_index = malloc( sizeof( PIOINT) * read_sample_count);
      assert( (data->pixel_index != NULL) && "Not enough memory to allocate pixel index array");
      assert( noDMC_readObject_PIOINT( Param->TOI_pixel_index, read_first_sample, read_sample_count, data->pixel_index) >= 0);
    }
    PIOINT *pixel_index = data->pixel_index;

    // prepare parameters for stim modules
    stimPar->BeginRing                  = loop_begin_ring; // first ring number to process
    stimPar->EndRing                    = loop_end_ring;   // last ring number to process
    stimPar->beginring_first_sample     = loop_first_sample - read_first_sample; // relative to buffer start
    stimPar->endring_last_sample        = loop_last_sample - read_first_sample;  // relative to buffer start
    stimPar->signal_first_sample_number = read_first_sample;
    stimPar->total_signal_length        = read_sample_count;


//==============================================================================
// systematics simulation

    // initialise noise TOI
    PIOFLOAT *noise = malloc( read_sample_count * sizeof( PIOFLOAT));
    assert( (noise != NULL) && "Not enough memory to allocate noise array");
    for (i = 0; i < read_sample_count; i++) {
      noise[i] = 0.0;
    }

    // add (photonic) one over f and white noise before convolution
    create_rngstreams( stimPar);
    if (Param->do_photonic_noise != 0.0) {
      assert( add_oof_and_white_noise( stimPar, noise, PHOTONIC_NOISE) == 0);
    }
    delete_rngstreams( stimPar);

    // bolometer time response (LFER) convolution
    if (Param->do_LFER_convolve == 1) {
      stimPar->conv_CONVOLVE = CONVOLVE;
      strcpy( stimPar->conv_R_FILTER,   "NONE");
      strcpy( stimPar->conv_LPF_FILTER, "NONE");
      assert( convolve( stimPar, noise) == 0);
    }
  
    // add (electronic) white noise after convolution
    create_rngstreams( stimPar);
    if (Param->do_electronic_noise != 0.0) {
      assert( add_oof_and_white_noise( stimPar, noise, ELECTRONIC_NOISE) == 0);
    }
    delete_rngstreams( stimPar);

    // tweak noise over 3 sigma before adding BL/OOF noise
    if (Param->do_despike_flag != 0.0) {
      assert( despike_flag( stimPar, noise) == 0);
    }

    // bolometer time response deconvolution and lowpass filtering
    if (Param->do_LFER_deconvolve == 1) {
      stimPar->conv_CONVOLVE = DECONVOLVE;
      if (brimo->freq >= 545) {
        strcpy( stimPar->conv_R_FILTER,   "COSINE");
        strcpy( stimPar->conv_LPF_FILTER, "NONE");
        // explicit 3pt filtering for 545 and 857 instead of lowpass filtering with convolve()
        assert( real_filter( stimPar, noise) == 0);
      } else {
        // freq < 545, low-pass filter inside tau_deconv
        strcpy( stimPar->conv_R_FILTER,   "NONE");
        strcpy( stimPar->conv_LPF_FILTER, "GAUSSCOSSQR");
      }
      if (strcmp( Param->deconv_LPF_FILTER, "0") != 0) {
        // use value from param file instead of default behavior
        strcpy( stimPar->conv_LPF_FILTER, Param->deconv_LPF_FILTER);
      }
      if (strcmp( Param->deconv_R_FILTER, "0") != 0) {
        // use value from param file instead of default behavior
        strcpy( stimPar->conv_R_FILTER, Param->deconv_R_FILTER);
      }
      // run deconvolution
      assert( convolve( stimPar, noise) == 0);
    }

    // add convolved/deconvolved/filtered noise to signal
    for (i = 0; i < read_sample_count; i++) {
      signal[i] += noise[i];
    }
    free( noise);

    // convert signal from KCMB to ADU, add baseline, add ADC nonlinearity and 4K lines, modulate
    create_rngstreams( stimPar);
    if (Param->do_sim_adu == 1) {
      assert( sim_adu( stimPar, signal) == 0);
    }
    delete_rngstreams( stimPar);

    // simulate compression / decompression step of transmission
    if (Param->do_compress_decompress == 1) {
      assert( compress_decompress( stimPar, signal) == 0);
//      assert( desire_codec( stimPar, signal) == 0);
    }

//==============================================================================
// systematics correction (aka TOI processing)

    // ADC nonlinearity correction and (some) 4K lines removal
    if (Param->do_correct_adc == 1) {
      assert( correct_adc( stimPar, signal) == 0);
    }

    // convert modulated signal from ADU to Volts and produce an extra 3pt demodulated output for despike
    if (Param->do_adu_to_volts == 1) {
      assert( adu_to_volts( stimPar, signal, NULL) == 0);
    }
  
    if (Param->do_bl_demod == 1) {
      assert( bl_demod( stimPar, signal) == 0);
    }
  
    // subtract glitch tails and thermal baseline from signal, convert from Volts to Watts  
    if (Param->do_gaindecorr == 1) {
      assert( gain_decorr( stimPar, signal, NULL) == 0);
    }
  
    // add oof_noise
    create_rngstreams( stimPar);
    if (Param->nharm_oof > 0) {
      assert( add_oof_noise( stimPar, signal) == 0);
    }
    delete_rngstreams( stimPar);


//==============================================================================
// optionally convert simulated TOI unit and combine it with an external TOI
// -  stimPar->calib contains KCMB2WATT, KCMB2WATT contains the conversion factor from MJy/sr to Watts for 545/857
// -  in 100-353 standard simulations, input TOIs are in KCMB and are converted to ADU, then Volts, then Watts
// -  in 545/857 simulations, input TOIs are in MJy/sr and not converted further (no ADU/ADC). So to have the
//    output TOIs in Watts, we must set do_calibrate = -1 to convert them from MJy/sr to Watts

    if (Param->do_calibrate == -1.0) {
      // convert signal from KCMB (MJy/sr) to Watts
      Param->do_calibrate = stimPar->calib;
    }
    if (Param->do_calibrate == 1.0) {
      // convert signal from Watts to KCMB (MJy/sr)
      Param->do_calibrate = 1.0 / stimPar->calib;
    }

    if (Param->add_final_toi_factor == 0.0) {
      strcpy( Param->add_final_toi, "0");
    }
    if ((Param->do_calibrate != 0) || (strcmp( Param->add_final_toi, "0") != 0)) {
      if (Param->do_calibrate == 0) {
        Param->do_calibrate = 1.0;
      } else {
        STIM_TRACE( 1, "multiplying final TOI by %g", Param->do_calibrate);
      }
      if (strcmp( Param->add_final_toi, "0") != 0) {
        if (Param->add_final_toi_factor == 0.0) {
          Param->add_final_toi_factor = 1.0;
        }
        if (data->add_toi == NULL) {
          data->add_toi = malloc( loop_sample_count * sizeof( PIOFLOAT));
          assert( data->add_toi != NULL);
          assert( noDMC_readObject_PIOFLOAT( Param->add_final_toi, loop_first_sample, loop_sample_count, data->add_toi) >= 0);
        }
        STIM_TRACE( 1, "adding TOI %g x %s", Param->add_final_toi_factor, Param->add_final_toi);
      }
      for (i = 0; i < loop_sample_count; i++) {
        signal[stimPar->beginring_first_sample + i] *= Param->do_calibrate;
        if (Param->add_final_toi_factor != 0.0) {
          signal[stimPar->beginring_first_sample + i] += Param->add_final_toi_factor * data->add_toi[i];
        }
      }
      if (!Param->stay_in_memory) {
        if (data->add_toi != NULL) {
          free( data->add_toi);
          data->add_toi = NULL;
        }
      }
    }


//==============================================================================
// projection of signal TOI in hpr_out


    if (Param->do_SHPR) {
      STIM_TRACE( 1, "projecting signal TOI to splined HPR");
      // splining HPRs of signals in watts (~1e-15) gives less numerical accuracy than
      // when signal is around 1.0. So when signal is in watts, we convert it to KCMB before splining it,
      // and back to Watts afterward

      for (int hpr_ring = loop_begin_ring; hpr_ring <= loop_end_ring; hpr_ring += 1) {       
        long beginringindex = BEGINRINGINDEX( hpr_ring);
        long endringindex   = ENDRINGINDEX(   hpr_ring);
        int ndata = endringindex - beginringindex + 1;
        
        long hpr_begin = (hpr_ring - begin_ring) * RINGSIZE;

        if (Param->do_calibrate == 0) {
          // convert to KCMB to have numerical values near 1.0
          for (i = 0; i < ndata; i++) {
            signal[beginringindex - read_first_sample + i] /= stimPar->calib;
          }
        }
  
        HPR_Cleaner( &pixel_index[beginringindex - read_first_sample], // in: one ring of {pixname}_HPRIDX_ABER_TotalFlag_dx11 TOI
                     &signal[beginringindex - read_first_sample],      // in: one ring of TOI to project to Splined HPR (SHPR)
                     &hpr_out[hpr_begin],                              // out: one produced SHPR ring of length RGSIZE, allocated by the caller
                     ndata,                                            // number of samples in <hidx> and <data>
                     NULL);
  
        if (Param->do_calibrate == 0) {
          // convert back to Watts
          for (i = 0; i < ndata; i++) {
            signal[beginringindex - read_first_sample + i] *= stimPar->calib;
          }
          for (i = 0; i < RINGSIZE; i++) {
            hpr_out[hpr_begin + i] *= stimPar->calib;
            assert( !isnan( hpr_out[hpr_begin + i]));
          }
        }
      }
    }
    else {
      STIM_TRACE( 1, "projecting signal TOI to legacy HPR");
      int pixel_hitcount[RINGSIZE];

      for (int hpr_ring = loop_begin_ring; hpr_ring <= loop_end_ring; hpr_ring += 1) {
        long beginringindex = BEGINRINGINDEX( hpr_ring);
        long endringindex   = ENDRINGINDEX(   hpr_ring);
        long hpr_begin      = (hpr_ring - begin_ring) * RINGSIZE;
  
        for (i = 0; i < RINGSIZE; i++) {
          hpr_out[hpr_begin + i] = 0.0;
          pixel_hitcount[i] = 0;
        }
  
        for (i = beginringindex - read_first_sample; i < endringindex + 1 - read_first_sample; i++) {
          if ((pixel_index[i] >= 0) && (pixel_index[i]) < RINGSIZE) {
            hpr_out[hpr_begin + pixel_index[i]] += signal[i];
            pixel_hitcount[pixel_index[i]] += 1;
          }
        }
  
        for (i = 0; i < RINGSIZE; i++) {
          if (pixel_hitcount[i] > 0) {
            hpr_out[hpr_begin + i] /= pixel_hitcount[i];
          }
          assert( !isnan( hpr_out[hpr_begin + i]));
        }
      }
    }

    // optionally write output signal and flag if Param->signal_out != "0"
    if (strcmp( Param->signal_out, "0") != 0) {
      sprintf( outname, "%s_%03d", Param->signal_out, (int)Param->random_seed);
      if (stimPar->mpi_rank == 0) {
        printf( "%s writing signal TOI to %s, samples (%ld-%ld)\n", stimPar->msg_prefix, Param->signal_out, loop_first_sample, loop_last_sample);
      }
      assert( BFT_writeObject( outname, "float32", loop_first_sample, loop_sample_count, signal + stimPar->beginring_first_sample) >= 0);
    }

    free( signal);
    if (!Param->stay_in_memory) {
      free( data->pixel_index);
      data->pixel_index = NULL;
    }
  } // end of loop on rings


//==============================================================================
// add a list of HPRs to stim HPR

  if (Param->flag_add_hpr_name == 1) {
    if (hpr_out == NULL) {
      if (stimPar->mpi_rank == 0) {
        STIM_TRACE(0, "stimexe WARNING: when hpr_out=0, add_hpr_* is ignored");
      }
    } else {
      if (data->add_hpr == NULL) {
        data->add_hpr = malloc( Param->n_add_hpr_name * sizeof( PIOFLOAT*));
        assert( data->add_hpr != NULL);
        for (int hpr_idx = 0; hpr_idx < Param->n_add_hpr_name; hpr_idx++) {
          data->add_hpr[hpr_idx] = malloc( rings_to_process * RINGSIZE * sizeof( PIOFLOAT));
          assert( data->add_hpr[hpr_idx] != NULL);
          assert( noDMC_readObject_PIOFLOAT( Param->add_hpr_name[hpr_idx], begin_ring * RINGSIZE, rings_to_process * RINGSIZE, data->add_hpr[hpr_idx]) >= 0);
        }
      }
      for (int hpr_idx = 0; hpr_idx < Param->n_add_hpr_name; hpr_idx++) {
        STIM_TRACE( 1, "adding HPR: %g x %s", Param->add_hpr_factor[hpr_idx], Param->add_hpr_name[hpr_idx]);
        for (i = 0; i < rings_to_process * RINGSIZE; i++) {
          hpr_out[i] += Param->add_hpr_factor[hpr_idx] * data->add_hpr[hpr_idx][i];
        }
      }
      if (!Param->stay_in_memory) {
        for (int hpr_idx = 0; hpr_idx < Param->n_add_hpr_name; hpr_idx++) {
          free( data->add_hpr[hpr_idx]);
        }
        free( data->add_hpr);
        data->add_hpr = NULL;
      }
    }
  }

  // write HPR to disk if requested
  if (strcmp( Param->hpr_out, "0") != 0) {
    sprintf( outname, "%s_%03d", Param->hpr_out, (int)Param->random_seed);
    STIM_TRACE( 1, "writing signal HPR to %s, samples (%ld-%ld)", outname, begin_ring * RINGSIZE, (begin_ring + rings_to_process) * RINGSIZE - 1);
    assert( BFT_writeObject( outname, "float32", begin_ring * RINGSIZE, rings_to_process * RINGSIZE, hpr_out) >= 0
            && "stim.c: error while writing <hpr_out> to disk.");
  }

  // trace rank walltime
  STIM_TRACE( 1, "%s rings %d-%d, walltime: %s", Param->bolometer, begin_ring, begin_ring + rings_to_process - 1, get_elapsed_time( &t0));

  return( 0);
}
