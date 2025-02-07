/*
 * Copyright 2018 Oticon A/S
 *
 * SPDX-License-Identifier: Apache-2.0
 */
/**
 * Module of the channel2G4 that takes care of the path loss calculations
 */

#include <ctype.h>
#include <math.h>
#include "bs_types.h"
#include "bs_tracing.h"
#include "bs_oswrap.h"
#include "channel_2G4Indoorv1_argparse.h"
#include "channel_2G4Indoorv1_int.h"

#define UNINITIALIZED 0
#define CONSTANT_ATT  1
#define FROM_FILE     2


#define MAXLINESIZE 2048

typedef struct { /*Structure to keep the status of each path file*/
  char * filename;
  FILE * fileptr;
  bs_time_t LastTime;
  bs_time_t NextTime;
  double LastDistValue;
  double LastPathLossValue;
  double NextDistValue;
} file_status_t;

typedef struct { /*All the status of the channel path loss*/
  uint ndevices;

  double distance_exp;

  double attenuation;
  double atxtra;
  char *matrix_file_name;

  char *path_mode; //For each of the NxN paths, in which mode does it run
  double *paths_att;

  file_status_t *distancefiles_status;
} pathloss_status_t;

static pathloss_status_t PL_status;

/**
 * Read a line from a file into a buffer (s), while
 * removing duplicate spaces (unless they are quoted), comments (#...), and ":",
 * and empty lines
 * The string will be null terminated
 *
 * Return: The number of characters copied into s
 */
static int PL_ReadLine(char *s, int size, FILE *stream){
  int c = 0, i=0;
  bool was_a_space = true;
  bool in_a_string = false;

  while ((i == 0) && (c != EOF)) {
    while ((i < size - 1) && ((c=getc(stream)) != EOF) && c!='\n') {
      if (isspace(c) && (!in_a_string)) {
        if (was_a_space) {
          continue;
        }
        was_a_space = true;
      } else {
        was_a_space = false;
      }
      if ((c == ':') && (!in_a_string)) {
        continue;
      }
      if (c=='"') {
        in_a_string = !in_a_string;
      }
      if (c == '#') {
        bs_skipline(stream);
        break;
      }
      s[i++] =c;
    }
  }
  s[i] = 0;

  if (i >= size - 1) {
    bs_trace_warning_line("Truncated line while reading from file after %i chars\n",size-1);
  }
  return i;
}

/*
 * Return the "length" of the string.
 * where the end is a "\0" or "\""
 */
static int PL_strlen(char *input){
  uint i=1;
  while (input[i] != 0 && input[i] != '"') {
    i++;
  }
  return i;
}

/**
 * Copy a string ("<string>") into output and null terminate it
 * note that the "" will be removed, and that the input pointer is expected to point at the first "
 */
static int PL_copy_string(char *input, char* output, uint *n){
  uint i=1;
  while (input[i] != 0 && input[i] != '"') {
    output[i-1] = input[i];
    i++;
  }
  *n = i;
  output[i-1] = 0;

  if (i == 1) {
    return -1;
  } else {
    return 0;
  }
}


/**
 * Calculate the average path loss in dB given a distance in meters
 */
static double PathLossFromDistance(double distance){
  double PL;
  if ( distance <= 0.0 ){
    bs_trace_warning_line("distance between devices = %f, this seems like an error..\n",distance);
  }

  PL = PL_status.distance_exp*10.0*log10(distance) + 39.60422483423212045872;
  //Ltotal = 20*log10f + N*log10d - 28
  //20*log10(2.4e3) - 28 = 39.60422483423212045872

  if ( PL < 20.0 ){
    static int never_complained_about_silly_distance = 1;
    if ( never_complained_about_silly_distance ){
      never_complained_about_silly_distance = 0;
      bs_trace_warning_line("distance between devices is very small (%.3fm).[the pathloos (%.1fdB) has been limited to 20dB] This channel model does not model near field conditions. Are you sure you want to do this? (this warning won't appear anymore)\n",distance, PL);
    }

    PL = 20;
  }

  return PL;
} //PathLossFromDistance()


/**
 * Initialize the reading/interpolation of path loss
 * from a distance file (for one path)
 */
static uint PL_Initiate_DistanceFile(file_status_t *this_status, uint index){
  //Try to read first line
  // if succeeded set LastTime, LastDistance to that, calculate LastPathLossValue
  //Try to read next line
  // if succeeded
  //    set NextTime and NextDistance to that
  //    return FROM_FILE;
  // if not succed (only 1 line in file)
  //    set the paths_att[index] = LastPathLossValue + atxtra
  //    return CONSTANT_ATT

  int read;
  char line_buf[MAXLINESIZE];
  char *filename = this_status->filename;
  this_status->fileptr = bs_fopen(filename, "r");

  //we try to read the first line:
  read = PL_ReadLine(line_buf, MAXLINESIZE, this_status->fileptr);
  if ( read == EOF ) {
    bs_trace_error_line("file %s seems empty\n",filename);
  }
  read = sscanf(line_buf, "%"SCNtime" %le", &this_status->LastTime , &this_status->LastDistValue);
  if (read < 2) {
    bs_trace_error_line("file %s seems corrupted\n",filename);
  } else {
    this_status->LastPathLossValue = PathLossFromDistance(this_status->LastDistValue);
  }

  //we try to read the second line:
  uint failed = 0;
  read = PL_ReadLine(line_buf, MAXLINESIZE, this_status->fileptr);
  if (read == 0) {
    failed = 1;
  } else {
    read = sscanf(line_buf, "%"SCNtime" %le", &this_status->NextTime , &this_status->NextDistValue);
    if (read == 0) {
      failed = 1;
    }
    if (read == 1) {
      bs_trace_error_line("file %s seems corrupted\n",filename);
    }
  }
  if (failed == 0) {
    return FROM_FILE;
  } else {
    PL_status.paths_att[index] = this_status->LastPathLossValue + PL_status.atxtra;
    return CONSTANT_ATT;
  }
}


/**
 * Return the path loss (for a path whose path loss is being read/interpolated from file data)
 */
static double PL_from_file(file_status_t *this_status, bs_time_t now, uint index){
  //while Now >= NextTime:
  //  move NextTime to LastTime
  //  move NextDistValue to LastDistValue, and calculate LastPathLossValue
  //  try to read next line
  //  if no next line
  //    overwrite the mode to CONSTANT_ATT, set the paths_att[index] = LastPathLossValue + atxtra
  //    return paths_att[index]
  //  else
  //    set NextTime and NextDistValue to the value,
  //
  //if Now <= LastTime:
  // return LastPathLossValue + PL_status.atxtra
  //elseif (in the middle)
  // interpolate distance,
  // calculate pathloss for that distance
  // return that pathloss + atxtra

  while (now >= this_status->NextTime) {
    uint failed = 0;
    int read;
    char line_buf[MAXLINESIZE];

    this_status->LastTime = this_status->NextTime;
    this_status->LastDistValue = this_status->NextDistValue;
    this_status->LastPathLossValue = PathLossFromDistance(this_status->LastDistValue);

    read = PL_ReadLine(line_buf, MAXLINESIZE, this_status->fileptr);
    if (read == 0) {
      failed = 1;
    } else {
      read = sscanf(line_buf, "%"SCNtime" %le", &this_status->NextTime , &this_status->NextDistValue);
      if (read == 0) {
        failed = 1;
      }
      if (read == 1) {
        bs_trace_error_line("file %s seems corrupted\n",this_status->filename);
      }
    }
    if ( failed == 1 ) { //no need to call this function again, the value wont change
      PL_status.paths_att[index] = this_status->LastPathLossValue + PL_status.atxtra;
      PL_status.path_mode[index] = CONSTANT_ATT;
      return PL_status.paths_att[index];
    }
  }


  if (now <= this_status->LastTime) {
    return this_status->LastPathLossValue + PL_status.atxtra;
  } else { //in between LastTime and NextTime
    double DistanceInterp = (this_status->NextDistValue - this_status->LastDistValue)
                            * ((double)(now-this_status->LastTime)
                            / (double)(this_status->NextTime - this_status->LastTime))
                            + this_status->LastDistValue;
    return PathLossFromDistance(DistanceInterp) + PL_status.atxtra;
  }

  //Unreachable
  return PL_status.atxtra;
}


static void ProcessMatrixFileLine(char *buf, uint *rx, uint *tx, double *distance, char **filename){
  uint read;
  uint off;
  uint ok;
  uint points_to_file = 0;
  char *endp;
  *filename = NULL;
  const char corruptmsg[] = "Corrupted distance matrix file (format txnbr rxnbr : {attenuation|\"filename\"}\n%s\n";

  *tx = strtoul(buf, &endp, 0);
  if (endp == buf) {
    bs_trace_error_line(corruptmsg, buf);
  }
  buf = endp;
  *rx = strtoul(buf, &endp, 0);
  if (endp == buf) {
    bs_trace_error_line(corruptmsg, buf);
  }
  buf = endp;

  { //check if next meaningfull char is " or a number
    off = 0;
    while (buf[off] != 0) {
      if (buf[off] == '"') {
        points_to_file = 1;
        break;
      }
      if ((buf[off] >= '0') && (buf[off] <= '9')) {
        break;
      }
      off++;
    } //while
  }
  if (points_to_file) {
    uint length = PL_strlen(&buf[off]);
    *filename = (char*) bs_calloc(length, sizeof(char));
    ok = PL_copy_string(&buf[off], *filename, &read);
    if (ok == -1) {
      bs_trace_error_line(corruptmsg, buf);
    }
  } else {
    read = sscanf(&buf[off], "%le", distance);
    if (read < 1) {
      bs_trace_error_line(corruptmsg, buf);
    }
  }
}


/*
 * Interfacing functions of the pathloss towards the rest of the channel:
 */


/*
 * Initialize all the path loss internals,
 * read the matrixfile (if it exists) and start loading the files it references to (if they exist)
 */
void InitPathLoss(ch_2G4I_args_t *args, uint n_devs){
  PL_status.ndevices    = n_devs;
  PL_status.attenuation = args->attenuation;
  PL_status.atxtra      = args->atxtra;
  PL_status.matrix_file_name = args->matrix_file_name;
  PL_status.distance_exp = args->distance_exp;

  PL_status.path_mode            = (char*)bs_calloc(n_devs*n_devs, sizeof(char));
  PL_status.paths_att            = (double*)bs_calloc(n_devs*n_devs, sizeof(double));
  PL_status.distancefiles_status = (file_status_t *)bs_calloc(n_devs*n_devs, sizeof(file_status_t));

  if ( args->matrix_file_name != NULL ){
    FILE *matrix_file;
    char line_buf[MAXLINESIZE];
    int read = 0;
    matrix_file = bs_fopen(args->matrix_file_name, "r");

    while (true) {
      uint rx, tx;
      double distance;
      char *filename;
      uint index;
      read = PL_ReadLine(line_buf, MAXLINESIZE, matrix_file);
      if ( read == 0 ) {
        break;
      }
      ProcessMatrixFileLine(line_buf, &rx, &tx, &distance, &filename);
      if ( ( rx >= n_devs ) || ( tx >= n_devs ) ){
        bs_trace_warning_line("The distances matrix file is trying to define the path from %i->%i, but only %i devices are set in the simulation => will be ignored\n",
                              tx, rx, n_devs);
        free(filename);
        continue;
      }
#if (AssumeChannelSimmetric == 1)
      if (rx >= tx) {
        bs_trace_warning_line("The distances matrix file is trying to define the path from %i->%i, but that one is assumed equal to %i->%i => will be ignored\n",
                              tx, rx, rx, tx);
        free(filename);
        continue;
      }
#endif

      index =  rx*PL_status.ndevices + tx;
      PL_status.distancefiles_status[index].filename = filename;

      if (filename == NULL){
        PL_status.paths_att[index] = PathLossFromDistance(distance) + PL_status.atxtra;
        if ( PL_status.path_mode[index] != UNINITIALIZED ) {
          bs_trace_warning_line("Redefinition of the path (%i->%i) path loss\n", tx, rx);
        }
        PL_status.path_mode[index] = CONSTANT_ATT;
      } else {
        PL_status.path_mode[index] = PL_Initiate_DistanceFile(&PL_status.distancefiles_status[index], index);

      }

    } //while
    fclose(matrix_file);
  }

  //set the default power level for all remaining paths
  uint tx,rx;
  for ( tx = 0 ; tx < n_devs; tx++){
#if (AssumeChannelSimmetric == 0)
    for ( rx = 0; rx < n_devs; rx++ ){
#else
    for ( rx = 0; rx < tx; rx++ ){
#endif
      uint index = rx*PL_status.ndevices + tx;
      if ( PL_status.path_mode[index] == UNINITIALIZED ) {
        if ( args->matrix_file_name != NULL ) {
          bs_trace_warning_line("The distance matrix file did not set the path %i->%i. It will be set to at + atxtra (%lf+%lf)\n",tx, rx, PL_status.attenuation, PL_status.atxtra);
        }
        PL_status.paths_att[index] = PL_status.attenuation + PL_status.atxtra;
        PL_status.path_mode[index] = CONSTANT_ATT;
      }
    }
  }
}


/**
 * Return the path loss from <tx> -> <rx> in this instant (<Now>)
 */
double CalculatePathLoss(uint tx, uint rx, bs_time_t now){
#if (AssumeChannelSimmetric == 1)
  if ( rx > tx ) {
    uint tmp = tx;
    tx = rx;
    rx = tmp;
  }
#endif
  uint index = rx*PL_status.ndevices + tx;
  if ( PL_status.path_mode[index] == CONSTANT_ATT ) {
    return PL_status.paths_att[index];
  } else if ( PL_status.path_mode[index] == FROM_FILE ) {
    return PL_from_file(&PL_status.distancefiles_status[index], now, index);
  } else {
    bs_trace_error_line("bad error\n");
    return 0;
  }
}


/**
 * Free all memory and close all files
 */
void clearPathLoss(){
  if ( PL_status.path_mode != NULL )
    free(PL_status.path_mode);
  if ( PL_status.paths_att != NULL )
    free(PL_status.paths_att);
  if ( PL_status.distancefiles_status != NULL ) {
    uint tx;
    uint rx;
    for (tx = 0 ; tx < PL_status.ndevices ; tx++){
      for ( rx = 0; rx < tx; rx++ ){
        uint index = rx*PL_status.ndevices + tx;
      if ( PL_status.distancefiles_status[index].fileptr != NULL ) {
        fclose(PL_status.distancefiles_status[index].fileptr);
      }
      if ( PL_status.distancefiles_status[index].filename != NULL ){
        free(PL_status.distancefiles_status[index].filename);
      }
      }
    }
    free(PL_status.distancefiles_status);
  }
}
