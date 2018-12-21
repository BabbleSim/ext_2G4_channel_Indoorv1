/*
 * Copyright 2018 Oticon A/S
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#include <math.h>
#include <complex.h>
#include <string.h>
#include <fftw3.h>
#include "bs_types.h"
#include "bs_oswrap.h"
#include "bs_utils.h"
#include "bs_tracing.h"
#include "bs_pc_2G4_types.h"
#include "bs_rand_main.h"
#include "bs_rand_inline.h"
#include "p2G4_pending_tx_rx_list.h"
#include "channel_2G4Indoorv1_argparse.h"
#include "channel_2G4Indoorv1_pathloss.h"
#include "channel_2G4Indoorv1_int.h"
#include "channel_if.h"

/*
 * This module provides a model of the 2G4 indoor channel (fading and path loss)
 *
 * Note that, although it internally calculates a plausible channel response (fading) of the whole ISM band over time,
 * it provides as output only the high level parameters used by the 2G4 Phy v1 as specified in channel_if.h
 *
 * The code consists mostly of 2 parts:
 *   * Path loss calculation part (see channel_2G4Indoorv1_pathloss.c & doc/README.path_loss.txt)
 *   * Fading calculation part. This is implemented in this file, see section 5 of
 *     the included document for documentation on the model
 *
 * The interfacing functions to the Phy are also implemented in this file
 */

//Constants for the fading calculation:

#define H_Ts (1.0/160e6)
///< channel impulse response sampling period
#define samplesPerMHz 2
///< 1 FFT tap per every 500kHz (we need to have a power of two number of samples per MHz => 1,2,4,8..
#define Log2SamplesPerMHz 1
#define Nsamples_ChannelAnalysis (160*samplesPerMHz)
///<How many taps for the whole 160 MHz BW (FFT of 160*2 taps in 160MHz)
#define Oversampling_recalc (16)
///<how often in time do we recalculate new taps compared to the Nyquist limit (for that given coherence time|DopplerSpread)

#define Doppler_filter_length  76
///Filter used to create the time correlation of the fading response
const static double Doppler_filter[] = {
   0.005656670561663,   0.002089101239094,   0.001275186438457,  -0.000617487188748,  -0.003646377236777 ,
  -0.007667389302124,  -0.012305735072649,  -0.016966079755007,  -0.020895584727681,  -0.023279764585737 ,
  -0.023360999807339,  -0.020607328895787,  -0.014873210087364,  -0.006433338933846,   0.003892512023787 ,
   0.014843174090853,   0.024808497608759,   0.032048601731793,   0.034956984169412,   0.032341891866916 ,
   0.023697155778454,   0.009405806866507,  -0.009178937572925,  -0.029773296857976,  -0.049344426393077 ,
  -0.064451784983745,  -0.071649222821988,  -0.067970674631071,  -0.051375931094812,  -0.021123173603960 ,
   0.021995350010165,   0.075614160898720,   0.135950397983927,   0.198150585855270,   0.256803759053764 ,
   0.306539769434757,   0.342645466746085,   0.361629024640802,   0.361629024640802,   0.342645466746085 ,
   0.306539769434757,   0.256803759053764,   0.198150585855270,   0.135950397983927,   0.075614160898720 ,
   0.021995350010165,  -0.021123173603960,  -0.051375931094812,  -0.067970674631071,  -0.071649222821988 ,
  -0.064451784983745,  -0.049344426393077,  -0.029773296857976,  -0.009178937572925,   0.009405806866507 ,
   0.023697155778454,   0.032341891866916,   0.034956984169412,   0.032048601731793,   0.024808497608759 ,
   0.014843174090853,   0.003892512023787,  -0.006433338933846,  -0.014873210087364,  -0.020607328895787 ,
  -0.023360999807339,  -0.023279764585737,  -0.020895584727681,  -0.016966079755007,  -0.012305735072649 ,
  -0.007667389302124,  -0.003646377236777,  -0.000617487188748,   0.001275186438457,   0.002089101239094 ,
   0.005656670561663}; //normalized to have a power gain of 0dB for a white gaussian input (=> has DC gain)
//TOOPT: this filter is quite an overkill and can be made smaller if needed for performance (the time taken by the channel is directly proportional to this filter length)

//Status of the fader:

static uint n_devices;
///<number of devices (N) => paths (NxN)
static ch_2G4I_args_t args;
///<command line arguments
static bs_time_t Recalc_time;
///< every how often do we need to recalculate the fader internal status
///<(we assume the fading process is practically static inside any Recalc_time or packet length)
static bs_time_t UncorrelatedTime;
///< after this amount of time the channel is fully uncorrelated (we can reset it instead of recalculating all intermediate steps)

static double rice_offset;
///<Ricean mean value of the first channel impulse response tap (h[0])
static uint N_h_taps;
///< number of channel impulse response taps == length(h[])
static double *H_taps_average_levels = NULL;
///<Average channel impulse response (tap level for each tap) for these given channel parameters (delay spread and rice k)
///< common for all paths
// [tap_number]
static double complex *H_Taps = NULL;
///<Current channel impulse response for this path
// [rx_nbr][tx_nbr][tap_number]
static double complex *ChannelResp = NULL;
///<Current channel frequency response (over the whole 160MHz)
// [rx_nbr][tx_nbr][Nsamples_ChannelAnalysis]

//Status of each tap filter:
static double complex *FiltersMemory = NULL;
///< Memory of the filters for the channel impulse response taps generation
///< For all NxN paths, and for each tap
///< the memory is kept here
// [rx_nbr][tx_nbr][Doppler_filter_tap][tap_number]
//Note the order of the last 2 indexes!! (Tap_number is the last index)
//we do it this way, so we iterate over taps when calculating to parallelize

static uint *filter_index = NULL;
///<pointer for the next input value into the filter memory. For each path (all taps of each path share this value)
//[rx][tx]
//End of fader very internal status

static bs_time_t *Last_Recalc_time = NULL;
///<Last time the H_taps and ChannelResp for this rx,tx has been recalculated
// [rx_nbr][tx_nbr]

/**********************************/
/* Wrapping of the fftw3 library: */
static double complex *FFTin, *FFTout;
static fftw_plan FFTplan;
static void PrepareFFT(){
  FFTin  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nsamples_ChannelAnalysis);
  memset(FFTin, 0, sizeof(fftw_complex)*Nsamples_ChannelAnalysis);
  FFTout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nsamples_ChannelAnalysis);
  FFTplan   = fftw_plan_dft_1d(Nsamples_ChannelAnalysis, FFTin, FFTout, FFTW_FORWARD, FFTW_ESTIMATE);
}
static void clearFFT(){
  fftw_destroy_plan(FFTplan);
  fftw_free(FFTin);
  fftw_free(FFTout);
  fftw_cleanup();
}
static inline void FFT_channel(double complex *in, double complex *out){
  memcpy(FFTin, in, N_h_taps * sizeof(complex double));
  //If N_h_taps would change on the fly we would need to clear the input buffer
  fftw_execute(FFTplan);
  memcpy(out, FFTout, Nsamples_ChannelAnalysis * sizeof(complex double));
}
/* End of wrapping of the fftw3 library */
/****************************************/

/**
 * Generate taps levels corresponding to a exponential delay spread
 * with rms delay spread = <RMS_DelaySpread>
 * and sampled at <Ts>
 *
 * <taps_levels> pointer to the array of Rayleigh taps levels [in natural units] (sigma)
 * they are normalized to have power 1 all together. [This corresponds to "Level tap_i" in the documentation]
 *
 * <N_taps> is the number of taps generated
 *
 * Note that taps_levels is allocated in this function
 */
static void Set_HTapsAverageLevels(uint *N_taps, double **taps_levels, double Ts, double RMS_DelaySpread){
  double MaxT = 3.5*RMS_DelaySpread;
  double time; //time in which we will calculate the taps
  uint tap;
  double rms_level = 0;

  *N_taps = BS_MAX(floor(MaxT / Ts),1);
  //Note that when copying to the input FFT buffer
  //we assume this value is constant

  *taps_levels = (double *)bs_calloc(*N_taps, sizeof(double) );
  for (tap = 0, time = 0; tap < *N_taps ; tap++ , time += Ts ){ //time is the time of each tap
    double tap_level = exp(-time/RMS_DelaySpread);
    (*taps_levels)[tap] = tap_level;
    rms_level += tap_level*tap_level;
  }

  rms_level = sqrt(rms_level);
  for (tap = 0; tap < *N_taps ; tap++ ){ //Normalize the taps to have power 1 all together
    (*taps_levels)[tap] /= rms_level;
  }

}

/**
 * Initialize the <n_devs>*<n_devs> faders internal status
 *
 * Note that all paths will have the same fading parameters, but each path will be an independent and uncorrelated process
 */
static void Channel_initFade(ch_2G4I_args_t *args, uint n_devs){
  //Initialize the faders
  uint tap;

  //we normalize the total Rice power to 1 = omega
  double taps_gain = sqrt( 1.0 / ( 2.0* ( 1.0 + args->Rice_K ) ) ); //== gausians std (not std^2) [In the documentation "R"]
  rice_offset = sqrt(args->Rice_K / (1.0 + args->Rice_K)); //in amplitude [In the documentation "nu (v)"]

  //Faders common part:
  double MaxDopplerSpread = args->DopplerSpeed/299.792458e6*2.4e9; //(approx 9Hz at 4km/h)

  if ( MaxDopplerSpread != 0 ) {
    Recalc_time = (1.0/MaxDopplerSpread)/(double)Oversampling_recalc*1e6; //Every how often do we need to recalculate the channel (Oversampling_recalc every coherence time)
  } else {
    Recalc_time = 1e6*1e6;
  }

  UncorrelatedTime = Recalc_time*Doppler_filter_length; //After this amount of time the channel will be fully independent of the previous state (not just "theoretically uncorrelated", as it is a full filter length farther)

  //taps level generator
  Set_HTapsAverageLevels(&N_h_taps, &H_taps_average_levels, H_Ts, args->RMS_DelaySpread);
  for (tap = 0 ; tap < N_h_taps; tap++){
    H_taps_average_levels[tap] *= taps_gain; //["tap_i sigma" in the documentation]
  }

  //Allocate arrays for NxN paths:
  FiltersMemory = (double complex*)bs_calloc( n_devs*n_devs*N_h_taps*Doppler_filter_length , sizeof(double complex) );
  filter_index = (uint *)bs_calloc( n_devs*n_devs, sizeof(uint) );

  H_Taps = (double complex*)bs_calloc( n_devs*n_devs*N_h_taps , sizeof(double complex));
  ChannelResp = (double complex*)bs_calloc( n_devs*n_devs*Nsamples_ChannelAnalysis , sizeof(double complex));

  Last_Recalc_time = (bs_time_t *)bs_calloc( n_devs*n_devs, sizeof(bs_time_t));

  {
    uint rx,tx;
    for (rx = 0; rx < n_devs; rx ++){
      for (tx = 0; tx < n_devs; tx ++){
        if (rx != tx){
          Last_Recalc_time[ rx*n_devs + tx ] = TIME_NEVER;
        } //if rx != tx
      } //for tx
    } //for rx
  }
  //Eof initialization of the faders

  PrepareFFT();
}



/**
 * Resets the fading internal status for the path tx -> rx
 * (internal status == filters' internal status)
 */
static inline void ResetPathInternalStatus(uint tx, uint rx){
  filter_index[ rx*n_devices + tx ] = 0;
  bs_random_Gaus_c_buffer(&FiltersMemory[ (rx*n_devices + tx)*Doppler_filter_length*N_h_taps ],
                      N_h_taps*Doppler_filter_length);
}

/**
 * Does one update of the fading internal status for the path tx -> rx
 * Equivalent to advancing that fader internal status by Recalc_time
 * Note that this only advances the filters' input/FIFO
 */
static inline void OneInternalFadeTick(uint tx, uint rx){
  uint filter_idx = filter_index[ rx*n_devices + tx ];

  //we load the new random number for this tap into this filter memory input index (filter_idx)
  bs_random_Gaus_c_buffer(&FiltersMemory[ ((rx*n_devices + tx)*Doppler_filter_length + filter_idx)*N_h_taps ],
                          N_h_taps);
  filter_idx = ( filter_idx + 1 ) % Doppler_filter_length;

  filter_index[rx*n_devices + tx] = filter_idx;
}

/**
 * Update the H_Taps for the path tx -> rx based on its fading internal status
 * (Note that it does not update the fading internal status, it only recalculates the taps)
 *
 * == Run N_h_taps filters in parallel (one for each tap)
 */
static inline void RecalculateHTaps(uint tx, uint rx){
  double complex *Taps = &H_Taps[ (rx*n_devices + tx)*N_h_taps ]; //we pick this path's Taps (filter's output)
  uint filter_idx = filter_index[ rx*n_devices + tx ]; //we pick this path's filters' status
  double complex acc[N_h_taps]; //One accumulator per filter
  int i;
  uint tap;

  memset(acc,0,N_h_taps*sizeof(double complex));

  for ( i = Doppler_filter_length -1 ; i >=0; i-- ){
    register double complex *FilterM = &FiltersMemory[ ((rx*n_devices + tx)*Doppler_filter_length + filter_idx)*N_h_taps ];
    register double df = Doppler_filter[i];
    for (tap = 0; tap < N_h_taps; tap++ ){
      acc[tap] += df*FilterM[tap];
    }
    filter_idx = ( filter_idx + 1 ) % Doppler_filter_length;
  }
  for (tap = 0 ; tap < N_h_taps; tap++){
    Taps[tap]  = acc[tap]*H_taps_average_levels[tap];
  }
  Taps[0] += rice_offset;

} //void ReCalculateMultipath()


/*
 * Recalculate the Frequency response of the ISM band (<ThisChResp>) given the time domain response (<ThisHTaps>)
 */
static inline void RecalculateChResp(double complex *ThisHTaps, double complex *ThisChResp){
  FFT_channel(ThisHTaps, ThisChResp);
}


/**
 * Advance (if necessary) the path <tx_i> -> <rx_i> until the time <Now>
 *
 * ( If necessary update the fading internal status for this path, and if so recalculate the frequency response
 *   If we are close enough in time to the last <This_Last_Recalc_time>, then no need to recalculate anything
 *  )
 */
static inline void UpdateFading(uint tx_i, uint rx_i, bs_time_t Now, bs_time_t *This_Last_Recalc_time, double complex *ThisChResp, double complex *ThisH_Taps){
  uint Recalced = 0;
  if ( (*This_Last_Recalc_time == TIME_NEVER) || ( Now - *This_Last_Recalc_time > UncorrelatedTime ) ) {
    ResetPathInternalStatus(tx_i, rx_i);
    *This_Last_Recalc_time = Now;
    Recalced = 1;
  } else {
    while ( *This_Last_Recalc_time + Recalc_time < Now ){
      OneInternalFadeTick(tx_i, rx_i);
      *This_Last_Recalc_time += Recalc_time;
      Recalced = 1;
    }
  }
  if (Recalced) {
    RecalculateHTaps(tx_i, rx_i); //this modifies ThisH_Taps
    RecalculateChResp(ThisH_Taps, ThisChResp);
  }
}




//Status kept (ONLY) between CalculateAveFade() and CalculateISISNR()
//This status is used only to save time to calculate the ISI, as CalculateISISNR can only be called right after CalculateAveFade()
// => CalculateISISNR() piggybacks on the internal calculations of CalculateAveFade()
static uint BW;   ///<BW in number of ChResp taps (channel resolution)
static int ChRespCenter;
static int FirstIndex, LastIndex;
static double D;

/**
 * Calculate the average fading for this tx -> rx path , for this modulation type and center frequency
 * Note that a positive fade is a gain, and a negative fade is an attenuation
 *
 * Where <ChResp> if the frequency response of the whole ISM band for this tx -> rx path
 *
 * Where ModulationType and CenterFreq are as defined in the Com IF
 * C_Modulation_t ModulationType: One of bs_COM_modulation*
 * C_Freq_t       CenterFreq    :From 0.0 to 80.0 in which frequency is the transmission centered
 *                     format 8.8 in MHz, offset relative to 2400 (can be negative for blockers < 2400 )
 */
static double CalculateAveFade(double complex *ChResp, p2G4_modulation_t ModulationType, p2G4_freq_t  CenterFreq){
  ///<This channel frequency response

  switch ( ModulationType ) {
    case P2G4_MOD_BLE:
      BW = 1*samplesPerMHz; break;
    case P2G4_MOD_PROP2M:
      BW = 2*samplesPerMHz; break;
    case P2G4_MOD_PROP3M:
      BW = 3*samplesPerMHz; break;
    case P2G4_MOD_PROP4M:
      BW = 4*samplesPerMHz; break;
    case P2G4_MOD_WLANINTER:
      BW = 16*samplesPerMHz; break;
    case P2G4_MOD_BLEINTER:
      BW = 1*samplesPerMHz; break;
    case P2G4_MOD_CWINTER:
      BW = 1; break;
    case P2G4_MOD_WHITENOISE1MHz:
      BW = 1*samplesPerMHz; break;
    case P2G4_MOD_WHITENOISE2MHz:
      BW = 2*samplesPerMHz; break;
    case P2G4_MOD_WHITENOISE4MHz:
      BW = 4*samplesPerMHz; break;
    case P2G4_MOD_WHITENOISE8MHz:
      BW = 8*samplesPerMHz; break;
    case P2G4_MOD_WHITENOISE16MHz:
      BW = 16*samplesPerMHz; break;
    case P2G4_MOD_WHITENOISE20MHz:
      BW = 20*samplesPerMHz; break;
    case P2G4_MOD_WHITENOISE40MHz:
      BW = 40*samplesPerMHz; break;
    case P2G4_MOD_WHITENOISE80MHz:
      BW = 80*samplesPerMHz; break;
    default:
      bs_trace_error_line("Unknown modulation type %u\n", ModulationType);
      break;
  }

  ChRespCenter = ( ( CenterFreq >> ( P2G4_freq_FRACB - Log2SamplesPerMHz - 1) ) + 1 ) >> 1; //Round CenterFreq to nearest tap of ChResp
  // P2G4_freq_FRACB is the number of fractional bits of CenterFreq
  // ( P2G4_freq_FRACB - Log2SamplesPerMHz ) : is the difference in number of fractional bits of the ChRespCenter and CenterFreq


  {
    /*
     Equivalent Matlab:
      ChannelResp = FFChannelResp( :, FirstFTap+1:LastFTap);
      N= LastFTap - FirstFTap;
      D = abs(sum(ChannelResp,2)/N).^2; %tap 0 power
    */
    int index;
    complex double acc = 0; // abs(sum(ChannelResp,2)/BW).^2; %tap 0 power
    FirstIndex = ChRespCenter - BW/2;
    LastIndex  =  ChRespCenter + BW/2 - 1; //from [-BW/2,BW/2) + Center
    LastIndex = BS_MAX(FirstIndex,LastIndex); //for the BW = 0 case

    if ( FirstIndex < -Nsamples_ChannelAnalysis ) {
      bs_trace_error_line("No transmitter can be lower than 2240MHz (really under 2360MHz is not a good idea already)\n");
    }
    if ( LastIndex >= Nsamples_ChannelAnalysis ) {
      bs_trace_error_line("No transmitter can be higher than 2560MHz (really over 2520MHz is not a good idea already)\n");
    }
    //The response of the channel is periodic every 160MHz, so we wrap it around and assume that the negative frequencies response is at 160MHz
    for ( index = FirstIndex ; (index < 0) && (index <= LastIndex) ; index ++ ){
      acc += ChResp[index + Nsamples_ChannelAnalysis];
    }
    for ( ; index <= LastIndex; index ++ ){
      acc += ChResp[index];
    }
    acc /= BW;

    D = creal(acc)*creal(acc) + cimag(acc)*cimag(acc); //abs(D)^2
  }

  return 10*log10(D); //AverageFadeLevel(:, ChannelCenter == ChannelCenters) = 10*log10(D);

} //double CalculateAveFade()


/**
 * Calculates the ISI level for a given channel and modulation scheme
 * given the Channel frequency response
 *
 * ** This function can only be called just after CalculateAveFade **
 * ** as it reuses its calculations and parameters (for execution speed sake)**
 */
static double CalculateISI(double complex *ChResp){
  /*
    Equivalent Matlab code:
      AllOther = sum(abs(ChannelResp).^2,2)/N - D; %all others taps powers
      ISI_SNR(:, ChannelCenter == ChannelCenters) = 10*log10(D./AllOther);
   */
  int index;
  double acc = 0;
  double AllOther;

  //The response of the channel is periodic every 160MHz, so we wrap it around and assume that the negative frequencies response is at 160MHz
  for ( index = FirstIndex ; (index < 0) && (index <= LastIndex) ; index ++ ){
    complex double a = ChResp[index + Nsamples_ChannelAnalysis];
    acc += creal(a)*creal(a) + cimag(a)*cimag(a); //sum(abs(ChannelResp).^2
  }
  for ( ; index <= LastIndex; index ++ ){
    complex double a = ChResp[index];
    acc += creal(a)*creal(a) + cimag(a)*cimag(a);
  }
  acc /= BW;

  AllOther = acc - D;

  if ( AllOther == 0 )
    return 100;
  else
    return 10*log10(D/AllOther);

}

/**
 * Clear Fading part internally allocated memory
 */
static void clearFading(){
  if ( H_taps_average_levels != NULL )
    free(H_taps_average_levels);
  if ( FiltersMemory != NULL )
    free(FiltersMemory);
  if ( filter_index != NULL )
    free(filter_index);
  if ( H_Taps != NULL )
    free(H_Taps);
  if ( ChannelResp != NULL )
    free(ChannelResp);
  if ( Last_Recalc_time != NULL )
    free(Last_Recalc_time);
  clearFFT();
}

/**
 * Channel IF functions:
 */

/**
 * Recalculate the fading and path loss of the channel in this current moment (<now>)
 * in between the N used paths and the receive path (<rxnbr>)
 *
 * inputs:
 *  tx_used    : array with n_devs elements, 0: that tx is not transmitting,
 *                                           1: that tx is transmitting,
 *               e.g. {0,1,1,0}: devices 1 and 2 are transmitting, device 0 and 3 are not.
 *  tx_list    : array with all transmissions status (the channel can check here the modulation type of the transmitter if necessary)
 *               (ignored in this channel)
 *  txnbr      : desired transmitter number (the channel will calculate the ISI only for the desired transmitter)
 *               (ignored in this channel)
 *  rxnbr      : device number which is receiving
 *               (ignored in this channel)
 *  now        : current time
 *               (ignored in this channel)
 *  att        : array with n_devs elements. The channel will overwrite the element i
 *               with the average attenuation from path i to rxnbr (in dBs)
 *               The caller allocates this array
 *  ISI_SNR    : The channel will return here an estimate of the SNR limit due to multipath
 *               caused ISI for the desired transmitter (in dBs)
 *               (This channel sets this value always to 100.0)
 *
 * Returns < 0 on error.
 * 0 otherwise
 */
int channel_calc(const uint *tx_used, tx_el_t *Tx_list, uint txnbr, uint rxnbr, bs_time_t Now, double *Atenuation, double *ISI_SNR){

  uint tx, tx_i;
  uint rx_i;
  for ( tx = 0 ; tx < n_devices; tx++ ){
    if ( tx_used[tx] ){
#if ( AssumeChannelSimmetric == 1 )
      if ( rxnbr > tx ) {
        tx_i = rxnbr;
        rx_i = tx;
      } else {
        rx_i = rxnbr;
        tx_i = tx;
      }
#else
      rx_i = rxnbr;
      tx_i = tx;
#endif

      if (tx != rxnbr){
        uint index = rx_i*n_devices + tx_i;
        bs_time_t *This_Last_Recalc_time = &Last_Recalc_time[index];
        double complex *ThisChResp = &ChannelResp[index*Nsamples_ChannelAnalysis];
        double complex *ThisH_Taps = &H_Taps[index*N_h_taps];

        UpdateFading(tx_i, rx_i, Now, This_Last_Recalc_time, ThisChResp, ThisH_Taps);

        Atenuation[tx] = CalculatePathLoss(tx_i,rx_i,Now) - CalculateAveFade(ThisChResp, Tx_list[tx].tx_s.radio_params.modulation, Tx_list[tx].tx_s.radio_params.center_freq);

        if ( tx == txnbr ){ //if this is the desired transmitter
          *ISI_SNR = CalculateISI(ThisChResp); //Calculate its ISI SNR
        }

      } else {
        bs_trace_error_line("A device transmitting while receiving is not supported\n");
      } //it tx != rx
    } //if tx_used
  } //for tx


  return 0;
}

/**
 * Initialize the channel
 */
int channel_init(int argc, char *argv[], uint n_devs){
  n_devices = n_devs;

  channel_2G4I_argparse(argc, argv, &args);
  InitPathLoss(&args, n_devices);
  Channel_initFade(&args, n_devs);

  return 0;
}

/**
 * Clean up: Free the memory the channel may have allocate
 * close any file descriptors etc.
 * (the simulation has ended)
 */
void channel_delete(){
  clearFading();
  clearPathLoss();
}


//#define TESTPL
#ifdef TESTFADE
// gcc -g -DTESTFADE *.c ../../../lib/libUtilv1.a ../../../lib/libRandv2.a -o channel_test -lfftw3 -lm -std=c99 -D_XOPEN_SOURCE=700 \
   -I ../../libUtilv1/src/ -I ../../libPhyComv1/src/ -I ../../ext_2G4_libPhyComv1/src/ -I ../../libRandv2/src/ -I ../../ext_2G4_phy_v1/src/
//./channel_test -args > test/params.m

/*
 * This test dumps to files the status of the fader so that it can be loaded in Matlab and checked
 *
 * This is a very hand made and run test, tweak in the test itself and check from the Matlab that you get what you expected
 */
int main(int argc, char**argv) {
  uint number_devices = 3;
  bs_time_t duration = 10e6;
  bs_time_t time = 0, cuttime = 50e6;
  FILE *CHRespF, *HF, *AveFfadeF;
  uint16_t Modulation = P2G4_MOD_PROP2M;
  uint cutted = 0;


  bs_random_init(23325);
  channel_init(argc-1, &argv[1], number_devices);
  printf("N_h_taps = %i ;\n",N_h_taps);
  printf("Recalc_time = %"PRItime" ;\n",Recalc_time);
  printf("Modulation = %u;\n",Modulation);
  printf("DelaySpread = %e;\n",args.RMS_DelaySpread);
  printf("DopplerSpeed = %e;\n",args.DopplerSpeed);
  printf("Rice_K = %e;\n",args.Rice_K);

  CHRespF = fopen("test/ChRespF.txt","w");
  HF = fopen("test/HF.txt","w");
  AveFfadeF = fopen("test/AveFade_ISI.txt","w");
  {
    uint tx_used[3] = {0,0,1};     uint txnbr = 2;     uint rxnbr = 1;
//    uint tx_used[3] = {0,1,0};     uint txnbr = 1;     uint rxnbr = 2;

    tx_el_t Tx_list[number_devices];
    Tx_list[2].tx_s.radio_params.center_freq = 0;
    Tx_list[2].tx_s.radio_params.modulation = P2G4_MOD_CWINTER;
    uint rx_i, tx_i;
    double Atenuation[number_devices];
    double ISI_SNR;
#if ( AssumeChannelSimmetric == 1 )
      if ( rxnbr > txnbr ) {
        tx_i = rxnbr;
        rx_i = txnbr;
      } else {
        rx_i = rxnbr;
        tx_i = txnbr;
      }
#else
      rx_i = rxnbr;
      tx_i = txnbr;
#endif

    while (time < duration){
      channel_calc(tx_used, Tx_list, txnbr, rxnbr, time, Atenuation, &ISI_SNR);
      time += Recalc_time;
      if ( cutted == 0 && time > cuttime ){
        time += UncorrelatedTime;
      }

      for (int i = 0; i < N_h_taps; i ++){
        double complex value = H_Taps[ (rx_i*n_devices + tx_i)*N_h_taps + i];
        fprintf(HF, "%lf %lf\n", creal(value), cimag(value));
      }
      for (int i = 0; i < Nsamples_ChannelAnalysis; i ++){
        double complex value =ChannelResp[(rx_i * number_devices + tx_i)*Nsamples_ChannelAnalysis + i];
        fprintf(CHRespF, "%lf %lf\n", creal(value), cimag(value));
      }


      uint CenterFreq = 0;
      double complex* ChResp = &ChannelResp[(rx_i * number_devices + tx_i)*Nsamples_ChannelAnalysis];
      for ( CenterFreq = 0  ; CenterFreq < 80 ; CenterFreq++ ) {
        double fade = CalculateAveFade(ChResp, Modulation, CenterFreq<<8);
        double ISI_SNR = CalculateISI(ChResp);

        fprintf(AveFfadeF,"%lf\t",fade);
        fprintf(AveFfadeF,"%lf\n",ISI_SNR);
      } //for

    }//while (time)
  }
  printf("BW = %i ;\n",BW/samplesPerMHz);

  fclose(AveFfadeF);
  fclose(CHRespF);
  fclose(HF);
  channel_delete();
}


#endif


#ifdef TESTPL
//gcc -g -DTESTPL *.c ../../../lib/libUtilv1.a ../../../lib/libRandv2.a -o channel_pl_test -lfftw3 -lm -std=c99 -D_XOPEN_SOURCE=700 \
 -I ../../libUtilv1/src/ -I ../../libPhyComv1/src/ -I ../../ext_2G4_libPhyComv1/src/ -I ../../libRandv2/src/ -I ../../ext_2G4_phy_v1/src/
//valgrind ./channel_pl_test -at=50  -atxtra=-10.1 #test with only constant
/*
valgrind  --leak-check=full ./channel_pl_test -at=50  -dist=test/trial_case_2.matrix #test with file with all cases
expected result:
time:      0      10      20      30      40      50      60      70      80      90     100     110
1->0    45.62   45.62   45.62   45.62   48.28   50.31   51.95   53.33   54.52   55.57   56.51   56.51
2->0    46.59   49.68   51.95   53.75   55.77   57.40   58.78   59.97   61.01   61.94   62.78   62.78
2->1    46.59   49.75   52.06   53.89   55.39   56.68   57.79   58.78   59.67   60.48   61.22   61.22
3->0    50.00   50.00   50.00   50.00   50.00   50.00   50.00   50.00   50.00   50.00   50.00   50.00
3->1    45.62   45.62   45.62   45.62   45.62   45.62   45.62   45.62   45.62   45.62   45.62   45.62
3->2    39.60   39.60   39.60   39.60   39.60   39.60   39.60   39.60   39.60   39.60   39.60   39.60
*/
/*
 * This test checks the path loss calculation part
 *
 * This is a very hand made and run test, tweak in the test itself and check from the Matlab that you get what you expected
 */
int main(int argc, char**argv) {
  uint number_devices = 4;
  bs_time_t duration = 120;
  bs_time_t time = 0;
  bs_time_t timeStep = 10;
  uint tx;
  uint rx;


  bs_random_init(23325);

  channel_2G4I_argparse(argc-1, &argv[1], &args);
  n_devices = number_devices;
  InitPathLoss(&args, n_devices);

  printf("time:\t");
  for ( time = 0 ; time < duration; time += timeStep) {
    printf("%4"PRItime"\t",time);
  }
  printf("\n");
  for (tx = 0 ; tx < 4; tx++){
    for ( rx = 0; rx<tx; rx++){
      printf("%u->%u\t",tx, rx);
      for ( time = 0 ; time < duration; time += timeStep) {
        double at =
            CalculatePathLoss(tx,rx,time);
        printf("%.2lf\t",at);
      }
      printf("\n");
    }
  }

  clearPathLoss();
}


#endif

