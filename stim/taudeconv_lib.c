
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "taudeconv_lib.h"

PIODOUBLE timeresp_pi = 3.141592653589793238462643383279502884197;
PIODOUBLE dummy=0;
PIOLONG write_file=0;

#define _PIOMALLOC(toto) malloc(toto)
#define _PIOFREE(toto) free(toto)


//  -------------------------------------------------------------------------
//
//  Name:  array_complex_mult
//  Purpose: Compute the complex multiplication of two arrays
//  Comments: 
//  Authors:  J.F. Macias-Perez
//  Revisions: 
//     7-mara-10: first easy version JFMP

// -------------------------------------------------------------------------
// Definition of variables 
void array_complex_mult(PIOLONG Ndata, PIODOUBLE *toi_1, PIODOUBLE *toi_2){
   PIOLONG n2,nj;
   PIODOUBLE *toistore;
   
  // allocation of  memory 
   n2 = (PIOLONG) (Ndata-1)/2; 
   toistore  = (PIODOUBLE *)_PIOMALLOC(Ndata*sizeof(PIODOUBLE));  

   for (long i=0; i < Ndata; i++) toistore[i]=toi_1[i]; 

    toi_1[0] = toistore[0]*toi_2[0]; 

    for (long j=1; j <= n2; j++){
   	 nj = Ndata-j;
     // complex multiplication (a,b) * (c,d) = (a*c - b*d , a*d+b*c)
     // real part is stored from 1:n2 ; imaginary part in reverse order from n2+2:Ndata-1
	 //       std::cout << j << toistore[j] << toi_2[j] << toistore[nj] << toi_2[nj] << std::endl;
         toi_1[j]  = toistore[j]*toi_2[j]  - toistore[nj]*toi_2[nj];
         toi_1[nj] = toistore[j]*toi_2[nj] + toistore[nj]*toi_2[j]; 
 
    }
    if (Ndata % 2 == 0) toi_1[n2+1] = toistore[n2+1]*toi_2[n2+1]; // imaginary part equal 0 for real data
    _PIOFREE(toistore);

} 

//  -------------------------------------------------------------------------
//
//  Name:  construct_timeresponse_inv
//  Purpose:  Give the inverse of the bolometer time response in Fourier space
//  Comments: 
//  Authors:  C. Renault, J. Aumont, J.F. Macias-Perez
//  Revisions: 
//     7-nov-05: first easy version JFMP
//     8-nov-05: checked order the fourier coefficients JFMP
//     9-nov-05: first working version
// -------------------------------------------------------------------------

void construct_timeresponse_inv(struct par_taudeconv pardeconv,PIOLONG Ndata,PIODOUBLE *tresponse){
   
   // define pi 

//  PIODOUBLE cpi2f2,a1a2,tau1tau2,a1tau2a2tau1;
  PIODOUBLE den;

  // getting the right method  
  PIODOUBLE *fr, *toi, *time;
//  PIODOUBLE normfilter;
  PIOLONG n21,n2;
  
  // constructing the frequency vector
  n21 = Ndata/2+1;
  n2=(Ndata-1)/2; 

  fr = (PIODOUBLE *)_PIOMALLOC(sizeof(PIODOUBLE)*n21);
  for (long i=0; i < n21; i++)  fr[i] = pardeconv.fsampling * (PIODOUBLE ) i / (PIODOUBLE ) Ndata;


  // Defining time vector

   toi  = (PIODOUBLE *)_PIOMALLOC(sizeof(PIODOUBLE)*Ndata);
   time  = (PIODOUBLE *)_PIOMALLOC(sizeof(PIODOUBLE)*Ndata);
 
   for (PIOLONG i=0; i < Ndata; i++) time[i] =( -1.0* (PIODOUBLE) Ndata + (PIODOUBLE ) i + 1.0) / pardeconv.fsampling;

/*
  if (!strcmp(pardeconv.method,"Time1cst")){ //if deconv method
      for (PIOLONG i=0 ; i < Ndata; i++) toi[i] =  pardeconv.A1 * exp(time[i] /  pardeconv.tau1 );
      normfilter = 0.0;
      for (PIOLONG i=0; i < Ndata; i++) normfilter+=toi[i];
      if (normfilter > 0.0) for (PIOLONG i=0; i < Ndata; i++) toi[i] = toi[i] * sqrt((PIODOUBLE ) Ndata) /normfilter;
      // -- Compute fft 
      c06eac(Ndata,toi, NAGERR_DEFAULT);
// -- Conjugate of the time constant
      for (long i=1; i <= n2 ; i++) toi[Ndata-i] *= -1.0;


  }else if (!strcmp(pardeconv.method,"Time2cst")){
 

      PIODOUBLE nmf1, nmf2;
      nmf1 = 0.0; nmf2=0.0; 
      for (PIOLONG i=0; i < Ndata; i++) {
	  nmf1 += exp(time[i] /  pardeconv.tau1 );
	  nmf2 += exp(time[i] /  pardeconv.tau2 );
      }
      for (PIOLONG i=0; i < Ndata; i++) toi[i] =  pardeconv.A1 * exp(time[i] /  pardeconv.tau1 )/nmf1  +
					           pardeconv.A2 * exp(time[i] /  pardeconv.tau2 )/nmf2;

      normfilter = 0.0;
      for (PIOLONG i=0; i < Ndata; i++) normfilter+=toi[i];
      if (normfilter > 0.0) for (PIOLONG i=0; i < Ndata; i++) toi[i] = toi[i] * sqrt((PIODOUBLE ) Ndata) /normfilter;
    
  
      // -- Compute fft 

      c06eac(Ndata,toi, NAGERR_DEFAULT);
// -- Conjugate of the time constant
      for (long i=1; i <= n2 ; i++) toi[Ndata-i] *= -1.0;

  }else if (!strcmp(pardeconv.method,"Time3cst")){

    //    tf_3time_constants(pardeconv,Ndata,time,toi);
    tf_3time_constants_fourier(pardeconv,Ndata,fr,toi);

// -----  Converting to standard Fourier transform convention
      for (PIOLONG index=1; index <= n2 ; index++) toi[Ndata-index] *= -1.0;


 }else if (!strcmp(pardeconv.method,"LFER4")){ 
 
// ------------- Compute  4 time constants TF 
     tf_4time_constants_fourier(pardeconv,Ndata,fr,toi); // Use JH old inverse convention
// -----  Converting to standard Fourier transform convention
     for (PIOLONG index=1; index <= n2 ; index++) toi[Ndata-index] *= -1.0;

// --- Analytical Model as computed by Jacques in the standard Fourier transform convention.
     tf_jh8(pardeconv,Ndata,fr, toi);  
     if (write_file) {
	 PIOSTRING filen = "/wrk/macias/tfJH8_direct_noclipping_v2.txt";
	 tf_write(filen, Ndata, fr, toi);
	 write_file = 0;
     }
     
// -- Check for zeros
     tf_modclipping(pardeconv,Ndata,toi);

 }else if (!strcmp(pardeconv.method,"LFER6")){ 
 
// ------------- Compute  6 time constants TF 
     tf_6time_constants_fourier(pardeconv,Ndata,fr,toi); // Use JH old inverse convention
// -----  Converting to standard Fourier transform convention
     for (PIOLONG index=1; index <= n2 ; index++) toi[Ndata-index] *= -1.0;

// --- Analytical Model as computed by Jacques in the standard Fourier transform convention.
     tf_jh8(pardeconv,Ndata,fr, toi);  
     if (write_file) {
	 PIOSTRING filen = "/wrk/macias/tfJH8_direct_noclipping_v2.txt";
	 tf_write(filen, Ndata, fr, toi);
	 write_file = 0;
     }
     
// -- Check for zeros
     tf_modclipping(pardeconv,Ndata,toi);

 }else

*/

 if (!strcmp(pardeconv.method,"LFER8")){ 
 
// ------------- Compute  8 time constants TF 
     tf_8time_constants_fourier(pardeconv,Ndata,fr,toi); // Use JH old inverse convention
// -----  Converting to standard Fourier transform convention
     for (PIOLONG index=1; index <= n2 ; index++) toi[Ndata-index] *= -1.0;

// --- Analytical Model as computed by Jacques in the standard Fourier transform convention.
     tf_jh8(pardeconv,Ndata,fr, toi);  
     if (write_file) {
	 PIOSTRING filen = "/wrk/macias/tfJH8_direct_noclipping_v2.txt";
	 tf_write(filen, Ndata, fr, toi);
	 write_file = 0;
     }
     
// -- Check for zeros
     tf_modclipping(pardeconv,Ndata,toi);
  
  }

// ----------------------------------------  
//    SHIFT IN TIME AND SPHASE SHIFT PART
// ----------------------------------------
    PIODOUBLE dumreal, dumimg, phase;
    for (long i=1; i <= n2 ; i++) {
    // complex multiplication (a,b) * (c,d) = (a*c - b*d , a*d+b*c)
    // ------ Using Jacques convention -------
      phase = -2.0 * timeresp_pi * fr[i] * pardeconv.global_offset / pardeconv.fsampling;
      dumreal = toi[i] * cos(phase) - toi[Ndata-i]*sin(phase);
      dumimg =  toi[i] * sin(phase) + toi[Ndata-i]*cos(phase);       
  
      toi[i]       = dumreal;
      toi[Ndata-i] = dumimg;
    }
    if (Ndata % 2 == 0) toi[n2+1] = toi[n2+1]; 
    //   std::cout << toi[n2+1] <<   std::endl;

//-----------------------------------------------------------
// ------ INVERTING THE FILTER FOR DECONVOLUTION ------------
//-----------------------------------------------------------
//	cout << "DOING CONVOLUTION ONLY !!!!" << endl;
    if ( pardeconv.convol){
	for (PIOLONG index=0; index < Ndata; index++) {
        if (isnan(toi[index])) {
          tresponse[index] = 0.0;
        }
        else {
          tresponse[index] = toi[index];
        }
      }
    }else{
//	cout << "DOING DECONVOLUTION ONLY !!!!" << endl;
	if (toi[0] != 0)  tresponse[0] = 1.0/toi[0];

	for (PIOLONG index=1; index <= n2 ; index++) {
	    den = (toi[index]*toi[index]+toi[Ndata-index]*toi[Ndata-index]);
	    tresponse[index] = toi[index]/den;
	    tresponse[Ndata-index] = -1.0* toi[Ndata-index]/den;

	}
	if (Ndata % 2 == 0) tresponse[n2+1] = 0.0; // WE SET DIRECTLY TO ZERO TO AVOID PROBLEMS!!

    }

  
//--   tf_write(Ndata, fr, tresponse);

// ---------  freeing out memory
  _PIOFREE(fr);
  _PIOFREE(toi);
  _PIOFREE(time);

}



//  -------------------------------------------------------------------------
//
//  Name:  construct_regularization_filter
//  Purpose:  
//  Comments: 
//  Authors:   J.F. Macias-Perez
//  Revisions: 
//     8-mars-10: first easy version JFMP
// -------------------------------------------------------------------------

void construct_lpf_filter(struct par_taudeconv pardeconv,PIOLONG Ndata,PIODOUBLE *tresponse){

  // getting the right method  
  PIODOUBLE *fr;
  PIOLONG n21,n2;
  
  // constructing the frequency vector
  n21 = Ndata/2+1;
  n2=(Ndata-1)/2; 

  fr = (PIODOUBLE *)_PIOMALLOC(sizeof(PIODOUBLE)*n21);
  for (long i=0; i < n21; i++)  fr[i] = pardeconv.fsampling * (PIODOUBLE ) i / (PIODOUBLE ) Ndata;
  
  if (!strcmp(pardeconv.lpf_method,"GAUSSCOSSQR")){ //if deconv method

//	tauparams->fgauss_lpf = Param->fgauss_lpf;
//	tauparams->fc_lpf = Param->fc_lpf;
//	tauparams->ffact_lpf = Param->ffact_lpf;
// fmax = fc + factor * (fsamp/2-fc)
     PIODOUBLE fmax =  pardeconv.fc_lpf + pardeconv.ffact_lpf * (pardeconv.fsampling/2.0-pardeconv.fc_lpf);
     PIODOUBLE xvar;
     tresponse[0] = 1.0;
     for (PIOLONG i=1; i <= n2 ; i++) {
	 tresponse[i] = 1.0;
	 tresponse[Ndata-i] = 0.0;
// Define a exp * cos^2 filter
// Bill and Brendan proposition
//
// We look at functions of the form: K(f) = Kgauss(f) Kcos (f)
//
// The Gaussian portion of the filter is Kgauss (f) = exp(-0.5 * (f/f_{gauss})^2)
//
//and the cosine portion of the filter is parameterized with $f_c$ and $factor$ as follows:
//
//    for f<fc: Kcos = 1
//    for f>fmax: Kcos = 0 where 
//    for f>fc>fmax: Kcos = cos^2 (pi/2 * x) where $x = (f-fc)/(fmax-fc)$ 

	 tresponse[i] = exp(-0.5 * fr[i] * fr[i]/pardeconv.fgauss_lpf/pardeconv.fgauss_lpf);

	 if (fr[i] < pardeconv.fc_lpf)  tresponse[i] = tresponse[i]* 1.0;
	 if (fr[i] > fmax)  tresponse[i] = 0.0;
	 if (  fr[i] > pardeconv.fc_lpf  && fr[i] < fmax ){
	     xvar = (fr[i] - pardeconv.fc_lpf)/(fmax - pardeconv.fc_lpf);
	     tresponse[i] = tresponse[i] *( 0.5 +  cos(  timeresp_pi* xvar)/2.0);
	 }
     }


     // JUST In CASE   
     if (Ndata % 2 == 0){
	  PIOLONG i = n2+1;

	 tresponse[i] = exp(-0.5 * fr[i] * fr[i]/pardeconv.fgauss_lpf/pardeconv.fgauss_lpf);

	 if (fr[i] < pardeconv.fc_lpf)  tresponse[i] = tresponse[i]* 1.0;
	 if (fr[i] > fmax)  tresponse[i] = 0.0;
	 if (  fr[i] > pardeconv.fc_lpf  && fr[i] < fmax ){
	     xvar = (fr[i] - pardeconv.fc_lpf)/(fmax - pardeconv.fc_lpf);
	     tresponse[i] = tresponse[i] *( 0.5 +  cos(  timeresp_pi* xvar)/2.0);
	 }

     }	  
  }

  // --- tf_write(Ndata, fr, tresponse);

  _PIOFREE(fr);
 

}

 
//  -------------------------------------------------------------------------
//
//  Name:  construct_regularization_filter
//  Purpose:  
//  Comments: 
//  Authors:   J.F. Macias-Perez
//  Revisions: 
//     8-mars-10: first easy version JFMP
// -------------------------------------------------------------------------

void construct_regularization_filter(struct par_taudeconv pardeconv,PIOLONG Ndata,PIODOUBLE *tresponse){

  // getting the right method  
  PIODOUBLE *fr;
  PIOLONG n21,n2;
  
  // constructing the frequency vector
  n21 = Ndata/2+1;
  n2=(Ndata-1)/2; 

  fr = (PIODOUBLE *)_PIOMALLOC(sizeof(PIODOUBLE)*n21);
  for (long i=0; i < n21; i++)  fr[i] = pardeconv.fsampling * (PIODOUBLE ) i / (PIODOUBLE ) Ndata;
  
  if (!strcmp(pardeconv.r_method,"COSINE")){ //if deconv method
      
     // define fcut filter assuming cosine from 0 at fcut to pi/2 at fnyquist
     PIODOUBLE fcut = fr[n2] - pardeconv.fwidth_rfilter;
     if (Ndata % 2 == 0)  fcut = fr[n2+1] - pardeconv.fwidth_rfilter;
     tresponse[0] = 1.0;
     for (PIOLONG i=1; i <= n2 ; i++) {
	 tresponse[i] = 1.0;
	 tresponse[Ndata-i] = 0.0;
// Define a cosine decay for the last part of the filter
// Xavier propose 1/2 * cos(x) + 1/2  where x = (f-fcut)/(fmax-fcut)  fmax-fcut = fwidth_rfilter
	 if (fr[i] >= fcut)  tresponse[i] = 0.5*cos(timeresp_pi*(fr[i]-fcut)/pardeconv.fwidth_rfilter)+0.5;
     }
     if (Ndata % 2 == 0) tresponse[n2+1] = 0.5*cos(timeresp_pi*(fr[n2+1]-fcut)/pardeconv.fwidth_rfilter)+0.5;
 
     // --- tf_write(Ndata, fr, tresponse);

// ------- GAUSSIAN OPTION ------
  }else if(!strcmp(pardeconv.r_method,"GAUSSIAN")){

// Define fcut filter assuming a Gaussian at 5 frwidth (sigma)
     PIODOUBLE fcut = fr[n2] - 5.0*pardeconv.fwidth_rfilter;
     if (Ndata % 2 == 0)  fcut = fr[n2+1] - 5.0* pardeconv.fwidth_rfilter;
     tresponse[0] = 1.0;
     for (PIOLONG i=1; i <= n2 ; i++) {
	 tresponse[i] =1.0;
	 tresponse[Ndata-i] = 0.0;
    // Define a gaussian decay for the last part of the filter
	 if (fr[i] >= fcut)  tresponse[i] = exp(-1.0*(fr[i]-fcut)*(fr[i]-fcut)/
					    2.0/pardeconv.fwidth_rfilter/pardeconv.fwidth_rfilter);

     }
    
     if (Ndata % 2 == 0){
	  PIOLONG i = n2+1;
	  tresponse[i] = exp(-1.0*(fr[i]-fcut)*(fr[i]-fcut)/
					    2.0/pardeconv.fwidth_rfilter/pardeconv.fwidth_rfilter);
     }	  
  }

  // --- tf_write(Ndata, fr, tresponse);

  _PIOFREE(fr);
 

}


// -------------------------------------------------------------------------//
void tf_8time_constants_fourier(struct par_taudeconv pardeconv, PIOLONG Ndata,PIODOUBLE *fr, PIODOUBLE *toi){
     PIODOUBLE normalization_frequency = 0.016; // TODO: make this a parameter
     PIODOUBLE tfnorm;
     PIODOUBLE zrimag,zrreal,omega;
     PIODOUBLE zrnorm1, zrnorm2 ,zrnorm3,zrnorm4,zrnorm5,zrnorm6,zrnorm7,zrnorm8;
     PIOLONG   n2=(Ndata-1)/2; 

// ------ COMPUTING NORMALIZATION FACTOR
     omega = 2.0*timeresp_pi*normalization_frequency;
     zrnorm1 = 1.0 + pardeconv.tau1*pardeconv.tau1*omega*omega;
     zrnorm2 = 1.0 + pardeconv.tau2*pardeconv.tau2*omega*omega;
     zrnorm3 = 1.0 + pardeconv.tau3*pardeconv.tau3*omega*omega;
     zrnorm4 = 1.0 + pardeconv.tau4*pardeconv.tau4*omega*omega;
     zrnorm5 = 1.0 + pardeconv.tau5*pardeconv.tau5*omega*omega;
     zrnorm6 = 1.0 + pardeconv.tau6*pardeconv.tau6*omega*omega;
     zrnorm7 = 1.0 + pardeconv.tau7*pardeconv.tau7*omega*omega;
     zrnorm8 = 1.0 + pardeconv.tau8*pardeconv.tau8*omega*omega;
     zrreal = pardeconv.A1/zrnorm1 +  pardeconv.A2/zrnorm2 +  pardeconv.A3/zrnorm3;
     if (pardeconv.A4 > 0.0) zrreal +=  pardeconv.A4/zrnorm4;
     if (pardeconv.A5 > 0.0) zrreal +=  pardeconv.A5/zrnorm5;
     if (pardeconv.A6 > 0.0) zrreal +=  pardeconv.A6/zrnorm6;
     if (pardeconv.A7 > 0.0) zrreal +=  pardeconv.A7/zrnorm7;
     if (pardeconv.A8 > 0.0) zrreal +=  pardeconv.A8/zrnorm8;
     zrimag = pardeconv.A1*pardeconv.tau1*omega/zrnorm1 + 
	 pardeconv.A2*pardeconv.tau2*omega/zrnorm2 + 
	 pardeconv.A3*pardeconv.tau3*omega/zrnorm3;
     if (pardeconv.A4 > 0.0 && pardeconv.tau4 > 0.0) zrimag +=  pardeconv.A4*pardeconv.tau4*omega/zrnorm4;
     if (pardeconv.A5 > 0.0 && pardeconv.tau5 > 0.0) zrimag +=  pardeconv.A5*pardeconv.tau5*omega/zrnorm5;
     if (pardeconv.A6 > 0.0 && pardeconv.tau6 > 0.0) zrimag +=  pardeconv.A6*pardeconv.tau6*omega/zrnorm6;
     if (pardeconv.A7 > 0.0 && pardeconv.tau7 > 0.0) zrimag +=  pardeconv.A7*pardeconv.tau7*omega/zrnorm7;
     if (pardeconv.A8 > 0.0 && pardeconv.tau8 > 0.0) zrimag +=  pardeconv.A8*pardeconv.tau8*omega/zrnorm8;
  
     tfnorm = sqrt( zrreal*zrreal + zrimag*zrimag);


  // ------ LOOP ON FREQUENCIES
    toi[0] = 1.0/tfnorm;
    for (long i=1; i <= n2 ; i++) {
	 omega = 2.0*timeresp_pi*fr[i];
	 zrnorm1 = 1.0 + pardeconv.tau1*pardeconv.tau1*omega*omega;
	 zrnorm2 = 1.0 + pardeconv.tau2*pardeconv.tau2*omega*omega;
	 zrnorm3 = 1.0 + pardeconv.tau3*pardeconv.tau3*omega*omega;
	 zrnorm4 = 1.0 + pardeconv.tau4*pardeconv.tau4*omega*omega;
	 zrnorm5 = 1.0 + pardeconv.tau5*pardeconv.tau5*omega*omega;
	 zrnorm6 = 1.0 + pardeconv.tau6*pardeconv.tau6*omega*omega;
	 zrnorm7 = 1.0 + pardeconv.tau7*pardeconv.tau7*omega*omega;
	 zrnorm8 = 1.0 + pardeconv.tau8*pardeconv.tau8*omega*omega;
	 zrreal = pardeconv.A1/zrnorm1 +  pardeconv.A2/zrnorm2 +  pardeconv.A3/zrnorm3;
	 if (pardeconv.A4 > 0.0) zrreal +=  pardeconv.A4/zrnorm4;
	 if (pardeconv.A5 > 0.0) zrreal +=  pardeconv.A5/zrnorm5;
	 if (pardeconv.A6 > 0.0) zrreal +=  pardeconv.A6/zrnorm6;
	 if (pardeconv.A7 > 0.0) zrreal +=  pardeconv.A7/zrnorm7;
	 if (pardeconv.A8 > 0.0) zrreal +=  pardeconv.A8/zrnorm8;
	 zrimag = pardeconv.A1*pardeconv.tau1*omega/zrnorm1 + 
	          pardeconv.A2*pardeconv.tau2*omega/zrnorm2 +  
	          pardeconv.A3*pardeconv.tau3*omega/zrnorm3;
	 if (pardeconv.A4 > 0.0 && pardeconv.tau4 > 0.0) zrimag += pardeconv.A4*pardeconv.tau4*omega/zrnorm4;
	 if (pardeconv.A5 > 0.0 && pardeconv.tau5 > 0.0) zrimag += pardeconv.A5*pardeconv.tau5*omega/zrnorm5;
	 if (pardeconv.A6 > 0.0 && pardeconv.tau6 > 0.0) zrimag += pardeconv.A6*pardeconv.tau6*omega/zrnorm6;
	 if (pardeconv.A7 > 0.0 && pardeconv.tau7 > 0.0) zrimag += pardeconv.A7*pardeconv.tau7*omega/zrnorm7;
	 if (pardeconv.A8 > 0.0 && pardeconv.tau8 > 0.0) zrimag += pardeconv.A8*pardeconv.tau8*omega/zrnorm8;
          
	 toi[i] = zrreal/tfnorm;
	 toi[Ndata-i] = zrimag/tfnorm;
     }

     if (Ndata % 2 == 0) {
	 PIOLONG i =n2+1;
	 omega = 2.0*timeresp_pi*fr[i];
	 zrnorm1 = 1.0 + pardeconv.tau1*pardeconv.tau1*omega*omega;
	 zrnorm2 = 1.0 + pardeconv.tau2*pardeconv.tau2*omega*omega;
	 zrnorm3 = 1.0 + pardeconv.tau3*pardeconv.tau3*omega*omega;
	 zrnorm4 = 1.0 + pardeconv.tau4*pardeconv.tau4*omega*omega;
	 zrnorm5 = 1.0 + pardeconv.tau5*pardeconv.tau5*omega*omega;
	 zrnorm6 = 1.0 + pardeconv.tau6*pardeconv.tau6*omega*omega;
	 zrnorm7 = 1.0 + pardeconv.tau7*pardeconv.tau7*omega*omega;
	 zrnorm8 = 1.0 + pardeconv.tau8*pardeconv.tau8*omega*omega;
	 zrreal = pardeconv.A1/zrnorm1 +  pardeconv.A2/zrnorm2 +  pardeconv.A3/zrnorm3;
         if (pardeconv.A4 > 0.0) zrreal +=  pardeconv.A4/zrnorm4;
         if (pardeconv.A5 > 0.0) zrreal +=  pardeconv.A5/zrnorm5;
         if (pardeconv.A6 > 0.0) zrreal +=  pardeconv.A6/zrnorm6;
         if (pardeconv.A7 > 0.0) zrreal +=  pardeconv.A7/zrnorm7;
         if (pardeconv.A8 > 0.0) zrreal +=  pardeconv.A8/zrnorm8;
	 zrimag = pardeconv.A1*pardeconv.tau1*omega/zrnorm1 + 
	          pardeconv.A2*pardeconv.tau2*omega/zrnorm2 + 
	          pardeconv.A3*pardeconv.tau3*omega/zrnorm3;
         if (pardeconv.A4 > 0.0 && pardeconv.tau4 > 0.0) zrimag +=  pardeconv.A4*pardeconv.tau4*omega/zrnorm4;
         if (pardeconv.A5 > 0.0 && pardeconv.tau5 > 0.0) zrimag +=  pardeconv.A5*pardeconv.tau5*omega/zrnorm5;
         if (pardeconv.A6 > 0.0 && pardeconv.tau6 > 0.0) zrimag +=  pardeconv.A6*pardeconv.tau6*omega/zrnorm6;
         if (pardeconv.A7 > 0.0 && pardeconv.tau7 > 0.0) zrimag +=  pardeconv.A7*pardeconv.tau7*omega/zrnorm7;
         if (pardeconv.A8 > 0.0 && pardeconv.tau8 > 0.0) zrimag +=  pardeconv.A8*pardeconv.tau8*omega/zrnorm8;
	 toi[n2+1]= zrreal/tfnorm;
     }    

     

}



//  -------------------------------------------------------------------------
//
//  Name:  tf_resonance_JH
//  Purpose:  Compute the resonance transfer function from JH
//  Comments: 
//  Authors:   J.F. Macias-Perez
//  Revisions: 
//     8-mars-10: first easy version JFMP
// -------------------------------------------------------------------------
void tf_resonance_JH(struct par_taudeconv pardeconv,PIOLONG Ndata, PIODOUBLE *fr, PIODOUBLE *toi){ 
     PIODOUBLE zrimag,zrreal,zrnorm, dreal,dimag, omega;
     PIOLONG n2=(Ndata-1)/2; 
     toi[0] = toi[0];
     for (long i=1; i <= n2 ; i++) {
	 omega = 2.0*timeresp_pi*fr[i];
	 zrnorm = (1.0 - pardeconv.res2*omega*omega)*(1.0 - pardeconv.res2*omega*omega)+ 
	           pardeconv.res3*omega*pardeconv.res3*omega;
	 zrreal = ((1.0 + pardeconv.res1*omega*omega)*(1-pardeconv.res2*omega*omega))/zrnorm;
	 zrimag = ((1.0 + pardeconv.res1*omega*omega)*pardeconv.res3*omega)/zrnorm;
	 dreal = toi[i]*zrreal - toi[Ndata-i]*zrimag;
	 dimag = toi[i]*zrimag + toi[Ndata-i]*zrreal;
	 toi[i] = dreal;
	 toi[Ndata-i] = dimag;
     }

     if (Ndata % 2 == 0) {
	 PIOLONG i =n2+1;
	 omega = 2.0*timeresp_pi*fr[i];
	 zrnorm = (1.0 - pardeconv.res2*omega*omega)*(1.0 - pardeconv.res2*omega*omega)+ 
	           pardeconv.res3*omega*pardeconv.res3*omega;
	 zrreal = ((1.0 + pardeconv.res1*omega*omega)*(1.0-pardeconv.res2*omega*omega))/zrnorm;
	 zrimag = ((1.0 + pardeconv.res1*omega*omega)*pardeconv.res3*omega)/zrnorm;
	 dreal = toi[i]*zrreal - toi[Ndata-i]*zrimag;
	 dimag = toi[i]*zrimag + toi[Ndata-i]*zrreal;

	 toi[n2+1]= dreal;
     }
}


//  -------------------------------------------------------------------------
//
//  Name:  tf_electronic_JH
//  Purpose:  Compute the electronic transfer function from JH
//  Comments: 
//  Authors:   J.F. Macias-Perez
//  Revisions: 
//     8-mars-10: first easy version JFMP
// -------------------------------------------------------------------------
void tf_electronic_JH(struct par_taudeconv pardeconv,PIOLONG Ndata, PIODOUBLE *fr, PIODOUBLE *toi){ 

    PIODOUBLE zeimag,zereal,zenorm, omega,dreal,dimag;
    PIOLONG n2=(Ndata-1)/2; 
    toi[0] = toi[0];
    for (long i=1; i <= n2 ; i++) {
	omega = 2.0*timeresp_pi*fr[i];
	zenorm = (1.0 - pardeconv.ele2 *omega*omega) * (1.0 - pardeconv.ele2 *omega*omega) + pardeconv.ele3*omega* pardeconv.ele3*omega;
	zereal = ((1.0- pardeconv.ele1*omega*omega)*(1.0-pardeconv.ele2*omega*omega))/zenorm;
	zeimag = ((1.0- pardeconv.ele1*omega*omega)*pardeconv.ele3*omega)/zenorm;
	dreal = toi[i]*zereal - toi[Ndata-i]*zeimag;
	dimag = toi[i]*zeimag + toi[Ndata-i]*zereal;
	toi[i] = dreal;
	toi[Ndata-i] = dimag;
      }

      if (Ndata % 2 == 0) {
	  PIOLONG i = n2+1;
	  omega = 2.0*timeresp_pi*fr[i];
	  zenorm = (1.0 - pardeconv.ele2 *omega*omega) * (1.0 - pardeconv.ele2 *omega*omega) + pardeconv.ele3*omega* pardeconv.ele3*omega;
	  zereal = ((1.0- pardeconv.ele1*omega*omega)*(1.0-pardeconv.ele2*omega*omega))/zenorm;
	  zeimag = ((1.0- pardeconv.ele1*omega*omega)*pardeconv.ele3*omega)/zenorm;
	  dreal = toi[i]*zereal - toi[Ndata-i]*zeimag;
	  dimag = toi[i]*zeimag + toi[Ndata-i]*zereal;

	  toi[n2+1]= dreal;
      }
}



// -------------------------------------------------------------------------
//
//  Name:  tf_elect_anal_jh8
//  Purpose:  Compute the tf with jh8
//  Comments: Based on JH's psuedo-code - this is how JH implements it, consistent with his definitions of parameters
//  Authors:   B.P.Crill
//  Revisions: 
//     8-juilliet-10: first easy version JFMP
// -------------------------------------------------------------------------
void tf_elect_anal_jh8(PIODOUBLE f, PIODOUBLE tau0, PIODOUBLE fmod,PIODOUBLE sphase,PIODOUBLE *tfr, PIODOUBLE *tfi){
 
  PIOLONG nn = 5; // order to which to exapand square wave bias
  
  PIODOUBLE omegamod,omegap,omegam;
  PIODOUBLE omega = 2.0 * timeresp_pi * f;
  PIODOUBLE fangmod = fmod * 2.0 * timeresp_pi ;
  PIODOUBLE zf1plur,zf1plui,zf1minr,zf1mini;
  PIODOUBLE zSKplur,zSKplui,zSKminr,zSKmini;
  PIODOUBLE zh4plur,zh4plui,zh4mini,zh4minr;
  PIODOUBLE zfelpr,zfelpi,zfelmr,zfelmi;
  PIOINT i1;
  PIOINT signe = -1;
  
  PIODOUBLE zden1r, zden1i,zden2,zden3;
  PIODOUBLE z3oplur,z3oplui,z3ominr,z3omini;
  PIODOUBLE arg;
  
  PIODOUBLE  tau1   = 1.e3*100.e-9;
  PIODOUBLE  tau3   = 10.e3*10.e-9;
  PIODOUBLE  tau4   = 51.e3*1.e-6;
  PIODOUBLE  zz3    = 510.e3;
  PIODOUBLE  zz4    = 1e-6;
  PIODOUBLE  zx1    = 18.7e3;
  PIODOUBLE  zx2    = 37.4e3;
  
  (*tfr) =  0.0;
  (*tfi) =  0.0;
  
  for ( i1 = 1; i1 <= nn; i1 ++) {
    signe = - signe;
    omegamod = (2.0*i1 - 1.0)*fangmod;
    omegap = omega + omegamod;
    omegam = -(omegamod - omega);
    
    zfelpr = 0.0; zfelpi=0.0; zfelmr=0.0; zfelmi=0.0;

    // Resonance transfer function
    lowpass(omegap,tau0,&zfelpr,&zfelpi);
    lowpass(omegam,tau0,&zfelmr,&zfelmi);
    //     std::cout << "RESONANCE FILTER  " <<  zfelpr  <<  " " <<  zfelpi <<  " " << zfelmr  <<  " " <<  zfelmi  <<  std::endl;

   // Electronic rejection filter
    zf1plur = 1.0;  zf1plui = 0.5 * omegap * tau1;
    complex_division(&zf1plur,&zf1plui,1.0,omegap*tau1);
    
    zf1minr = 1.0;  zf1mini = 0.5 * omegam * tau1;
    complex_division(&zf1minr,&zf1mini,1.0,omegam*tau1);
    
    complex_multiplication(&zfelpr,&zfelpi,zf1plur,zf1plui);
    complex_multiplication(&zfelmr,&zfelmi,zf1minr,zf1mini);
    // std::cout << "REJECTION FILTER  " <<  zfelpr  <<  " " <<  zfelpi <<  " " << zfelmr  <<  " " <<  zfelmi  <<  std::endl;
    
    // Sallen-Key high pass filter
    zSKplur = 0.0; zSKplui = tau4*omegap;
    complex_division(&zSKplur,&zSKplui,1.0,omegap*tau4);
    zSKminr = 0.0; zSKmini = tau4*omegam;
    complex_division(&zSKminr,&zSKmini,1.0,omegam*tau4);
    
    complex_multiplication(&zfelpr,&zfelpi,zSKplur,zSKplui);
    complex_multiplication(&zfelmr,&zfelmi,zSKminr,zSKmini);
    
    // Square Sallen-Key term
    complex_multiplication(&zfelpr,&zfelpi,zSKplur,zSKplui);
    complex_multiplication(&zfelmr,&zfelmi,zSKminr,zSKmini);    
 
    //  std::cout << "SALLEN-KEY FILTER  " <<  zfelpr  <<  " " <<  zfelpi <<  " " << zfelmr  <<  " " <<  zfelmi  <<  std::endl;
  
    // Sign Reverse & gain
    complex_multiplication(&zfelpr,&zfelpi,-5.1,0.0);
    complex_multiplication(&zfelmr,&zfelmi,-5.1,0.0);
    
    // Single pole lowpass filter
    lowpass(omegap,tau3,&zh4plur,&zh4plui);
    zh4plur*=1.5;zh4plui*=1.5;
    lowpass(omegam,tau3,&zh4minr,&zh4mini);
    zh4minr*=1.5;zh4mini*=1.5;
    
    complex_multiplication(&zfelpr,&zfelpi,zh4plur,zh4plui);
    complex_multiplication(&zfelmr,&zfelmi,zh4minr,zh4mini);
    
//  std::cout << "Single pole FILTER  " <<  zfelpr  <<  " " <<  zfelpi <<  " " << zfelmr  <<  " " <<  zfelmi  <<  std::endl;
 //     // Third order equation


// plus frequency


//     zden3 = omegap*omegap*omegap*zx1*zx1*zz3*zx2*zx2*1.e-16*zz4;

//     // Corrected from the extra zx1*zx1*zx2*zx2*1.e-16 term !!!!
//     zden2 = omegap*omegap*(zx1*zx2*zx2*zz3*1.e-16 +
// 			   zx1*zx1*zx2*zx2*1.e-16+  
// 			   zx1*zx2*zx2*zz3*zz4*1.e-8);
    
//     // VERY WEIRD HERE
//     zden1 = omegap * (zx1*zx2*zx2*1.e-8+zx2*zz3*zx1*zz4) + zx2*zx1;
//     z3oplur = 0.0; z3oplui = 2.0*zx2*zx1*zz3*zz4*omegap;
//     complex_division(&z3oplur,&z3oplui,-1.0*zden2 + zx2*zx1 ,zden1-zden3);

   zden3 = -1.0 * omegap*omegap*omegap*zx1*zx1*zz3*zx2*zx2*1.0e-16*zz4;    
   zden2 = -1.0 * omegap*omegap*(zx1*zx2*zx2*zz3*1.e-16 +
 			   zx1*zx1*zx2*zx2*1.e-16+  
 			   zx1*zx2*zx2*zz3*zz4*1.e-8);
 

    
   zden1i = omegap * (zx1*zx2*zx2*1.e-8+zx2*zz3*zx1*zz4) + zden3; 
   zden1r = zx2*zx1 + zden2;
   
   // std::cout <<  "DENOMINATOR  PLUS" <<  zden1r << " " << zden1i << std::endl;
 
   z3oplur = 0.0; z3oplui= 2.0*zx2*zx1*zz3*zz4*omegap;
   complex_division(&z3oplur,&z3oplui, zden1r , zden1i);

   
    
// minus frequency
//     zden3 = omegam*omegam*omegam*zx1*zx1*zz3*zx2*zx2*1.e-16*zz4;


//     zden2 = omegam*omegam*(zx1*zx2*zx2*zz3*1.e-16 +
//        zx1*zx1*zx2*zx2*1.e-16+
//        zx1*zx2*zx2*zz3*zz4*1.e-8);
    
//     zden1 = omegam * (zx1*zx2*zx2*1.e-8+zx2*zz3*zx1*zz4) + zx2*zx1;
    
//     z3ominr = 0.0; z3omini = 2.0*zx2*zx1*zz3*zz4*omegam;
//     complex_division(&z3ominr,&z3omini,-1.0*zden2 + zx2*zx1,zden1-zden3);
    
   zden3 = -1.0 * omegam*omegam*omegam*zx1*zx1*zz3*zx2*zx2*1.0e-16*zz4;    
   zden2 = -1.0 * omegam*omegam*(zx1*zx2*zx2*zz3*1.e-16 +
 			         zx1*zx1*zx2*zx2*1.e-16+  
 			         zx1*zx2*zx2*zz3*zz4*1.e-8);
   zden1i = omegam * (zx1*zx2*zx2*1.e-8+zx2*zz3*zx1*zz4) + zden3; 
   zden1r = zx2*zx1 + zden2;
   
   // std::cout <<  "DENOMINATOR  MINUS" <<  zden1r << " " << zden1i << std::endl;

   z3ominr = 0.0; z3omini= 2.0*zx2*zx1*zz3*zz4*omegam;
   complex_division(&z3ominr,&z3omini, zden1r , zden1i);

    

// update the main variablex
   complex_multiplication(&zfelpr,&zfelpi,z3oplur,z3oplui);  
   complex_multiplication(&zfelmr,&zfelmi,z3ominr,z3omini);    
   //  std::cout << "BIG FILTER  " <<  zfelpr  <<  " " <<  zfelpi <<  " " << zfelmr  <<  " " <<  zfelmi  <<  std::endl;
    
 // ---- Averaging effect
    arg = timeresp_pi * omegap / (2.0*fangmod);
    complex_multiplication(&zfelpr,&zfelpi,-1.0*sin(arg)/arg,0.0);
    //  std::cout << "arg  " <<  arg  <<  "  " << -1.0*sin(arg)/arg << std::endl;
    arg = timeresp_pi * omegam / (2.0*fangmod);
    complex_multiplication(&zfelmr,&zfelmi,-1.0*sin(arg)/arg,0.0);
    // std::cout << "arg  " <<  arg  <<  "  " << -1.0*sin(arg)/arg << std::endl;
 
    complex_multiplication(&zfelpr,&zfelpi,cos(sphase*omegap), sin(sphase*omegap));
    complex_multiplication(&zfelmr,&zfelmi,cos(sphase*omegam), sin(sphase*omegam));
    //  std::cout << "COSP  " << cos(sphase*omegap)  <<  "  " << sin(sphase*omegap) << std::endl;
    // std::cout << "COSM  " << cos(sphase*omegam)  <<  "  " << sin(sphase*omegam) << std::endl;
    //  std::cout << "SPHASE  " <<  sphase  <<  std::endl;
    //   std::cout << "AVERAGE FILTER  " <<  zfelpr  <<  " " <<  zfelpi <<  " " << zfelmr  <<  " " <<  zfelmi  <<  std::endl;
 
  

    // Add to total
    (*tfr) +=  (signe/(2.0*i1-1))*(zfelpr+zfelmr);
    (*tfi) +=  (signe/(2.0*i1-1))*(zfelpi+zfelmi);
    //std::cout << "NORM  " <<  (*tfr)  <<  " " <<  (*tfi) <<  " "  <<  std::endl;
 

  }
  
  
  
}


//  -------------------------------------------------------------------------
//
//  Name:  lowpass
//  Purpose:  Compute lowpass filter
//  Comments: Based on Brendan code.
//  Authors:   J.F. Macias-Perez
//  Revisions: 
//     8-mars-10: first easy version JFMP
//  -------------------------------------------------------------------------

void lowpass(PIODOUBLE w, PIODOUBLE tau, PIODOUBLE *zr, PIODOUBLE *zi) {
  PIODOUBLE a = w * tau; 
  PIODOUBLE b = 1.0 + a * a;
  (*zr) = 1.0 / b ;
  (*zi) = - a / b ;
}


//  -------------------------------------------------------------------------
//
//  Name:  set_complex_exp
//  Purpose: set complex exponential 
//  Comments: Based on Brendan code.
//  Authors:   J.F. Macias-Perez
//  Revisions: 
//     8-mars-10: first easy version JFMP
//  -------------------------------------------------------------------------
void set_complex_exp(PIODOUBLE mod, PIODOUBLE arg,PIODOUBLE *zr,PIODOUBLE *zi) {
  (*zr) = mod * cos(arg);
  (*zi) = mod * sin(arg);
}


//  -------------------------------------------------------------------------
//
//  Name: complex_multiplication 
//  Purpose: complex multiplication
//  Comments: 
//  Authors:   J.F. Macias-Perez
//  Revisions: 
//     8-mars-10: first easy version JFMP
//  -------------------------------------------------------------------------
void complex_multiplication(PIODOUBLE *r1, PIODOUBLE *i1,PIODOUBLE r2,PIODOUBLE i2){
 
    PIODOUBLE dumr, dumi;
    
    dumr = (*r1) * r2 - (*i1) * i2;
    dumi = (*r1) * i2 + (*i1) * r2;
    (*r1) = dumr;
    (*i1) = dumi;
}

//  -------------------------------------------------------------------------
//
//  Name: complex_division
//  Purpose: complex division
//  Comments: 
//  Authors:   J.F. Macias-Perez
//  Revisions: 
//     8-mars-10: first easy version JFMP
//  -------------------------------------------------------------------------
void complex_division(PIODOUBLE *r1, PIODOUBLE *i1,PIODOUBLE r2,PIODOUBLE i2){
 
    PIODOUBLE dumr, dumi, den;
    den = r2 * r2 + i2 * i2;
    dumr = (*r1) * r2 + (*i1) * i2;
    dumi = -(*r1) * i2 + (*i1) * r2;
    (*r1) = dumr/den;
    (*i1) = dumi/den;
}



//  -------------------------------------------------------------------------
//
//  Name:  tf_jh8
//  Purpose:  Compute the resonance transfer function from JH
//  Comments: 
//  Authors:   B.P. Crill, J.F. Macias-Perez
//  Revisions: 
//     8-july-10: first easy version 
// -------------------------------------------------------------------------
void tf_jh8(struct par_taudeconv pardeconv,PIOLONG Ndata, PIODOUBLE *fr, PIODOUBLE *toi){ 
     PIODOUBLE tfr,tfi,tfr0,tfi0;
     PIODOUBLE dreal,dimag;
     PIOLONG n2=(Ndata-1)/2; 
     PIODOUBLE normalization_frequency = 0.016; // TODO: make this a parameter
     //PIODOUBLE normalization_frequency = 0.0; // to make compatible with Jacques function
     
     tfr0=0.0; tfi0=0.0;
     tf_elect_anal_jh8(normalization_frequency,pardeconv.tau_stray,pardeconv.fsampling/2.0, pardeconv.sphase,&tfr0,&tfi0);
 
//     std::cout << pardeconv.tau_stray  << " " << pardeconv.fsampling/2.0  << " " << pardeconv.sphase  << std::endl; 
//     std::cout << normalization_frequency << " " << tfr0  << " " << tfi0 << std::endl; 
      
     toi[0] = toi[0];  // Include normalization
     for (long i=1; i <= n2 ; i++) {
       tf_elect_anal_jh8(fr[i],pardeconv.tau_stray,pardeconv.fsampling/2.0, pardeconv.sphase,&tfr,&tfi);
       
    
   
// Do we really want to normalize this way    
       complex_division(&tfr,&tfi,tfr0,tfi0);  // We normalize the response to the zero frequency,.
       
       dreal = toi[i]*tfr - toi[Ndata-i]*tfi;
       dimag = toi[i]*tfi + toi[Ndata-i]*tfr;
       
      
       toi[i] = dreal;
       toi[Ndata-i] = dimag;
       //logger <<i << " " << fr[i]<<" "<< toi[i]<< PIO::endl;     
     }

     if (Ndata % 2 == 0){
	 tf_elect_anal_jh8(fr[n2+1],pardeconv.tau_stray,pardeconv.fsampling/2.0, pardeconv.sphase,&tfr,&tfi);	 
         toi[n2+1]= tfr/tfr0;
     }
     
}


//  -------------------------------------------------------------------------
//
//  Name:  tf_modclipping
//  Purpose:  Compute the resonance transfer function from JH
//  Comments: 
//  Authors:   J.F. Macias-Perez
//  Revisions: 
//     8-mars-10: first easy version JFMP
// -------------------------------------------------------------------------

void tf_modclipping(struct par_taudeconv pardeconv, PIOLONG Ndata, PIODOUBLE *toi){

      PIODOUBLE checknorm;
      PIOLONG n2 = (Ndata-1)/2; 
     
      for (PIOLONG i=1; i <= n2 ; i++) {
	  checknorm = sqrt(toi[i]*toi[i]+toi[Ndata-i]*toi[Ndata-i]);
	  //        if (i > n2-8) std::cout << "NORM " << checknorm << std::endl;
	  if ((checknorm < pardeconv.tfmodclip) & (checknorm != 0.0)){
	      toi[i] =  toi[i]/checknorm* pardeconv.tfmodclip;
	      toi[Ndata-i] =  toi[Ndata-i]/checknorm* pardeconv.tfmodclip;
	      //             if (i > n2-8) std::cout << "CLIPIN " << toi[i]  << " " << toi[Ndata-i]  << std::endl;

	  }else if (checknorm == 0.0) {
	      toi[i] =  pardeconv.tfmodclip;
	      toi[Ndata-i] = 0.0;
	  }
      }

       if (Ndata % 2 == 0){
              
	    if (sqrt(toi[n2+1]*toi[n2+1]) < pardeconv.tfmodclip)   toi[n2+1] =  pardeconv.tfmodclip;
       }
}




//  -------------------------------------------------------------------------
//
//  Name:  tf_write
//  Purpose:  write tf for testing purposes
//  Comments: 
//  Authors:   J.F. Macias-Perez
//  Revisions: 
//     8-mars-10: first easy version JFMP
// -------------------------------------------------------------------------

void tf_write( PIOSTRING filename, PIOLONG Ndata, PIODOUBLE *fr, PIODOUBLE *toi){
/*
    if (dummy == 0) {
      PIOLONG n21,n2;
      n21 = Ndata/2+1;
      n2 = (Ndata-1)/2; 
      

      //     FILE* f=fopen("/wrk/macias/tfLFER10.txt", "r");  
     
     ofstream myfile; 
     //    myfile.open ("/wrk/macias/tfJH10_direct_noclipping.txt");
     myfile.open (filename);
     myfile <<  fr[0]  << " "  <<  toi[0] << " " << 0.0 << "\n";
//       fprintf(f, "%g \t  %g \t  %g \n", fr[0],toi[0],0.0);  
     for (PIOLONG i=1; i <= n2 ; i++) {
	  myfile <<  fr[i]  << " " <<  toi[i] << " " << toi[Ndata-i] << "\n";
//          fprintf(f, "%g \t  %g \t  %g \n", fr[i],toi[i],toi[Ndata-i]);    // Print "x=3"  Other conversions:

       }

       if (Ndata % 2 == 0){
        myfile <<  fr[n2+1]  << " " <<  toi[n2+1] << " " << 0.0 << "\n";
//          fprintf(f, "%g \t  %g \t  %g \n", fr[n2+1],toi[n2+1],0.0);    // Print "x=3"  Other conversi           
       }
       myfile.close();
//       fclose(f);                // Close file f
       dummy =1;
  
    }
*/
}

