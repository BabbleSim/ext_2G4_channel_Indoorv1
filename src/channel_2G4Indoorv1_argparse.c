/*
 * Copyright 2018 Oticon A/S
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdint.h>
#include <inttypes.h>
#include <limits.h>
#include <stddef.h>
#include "channel_2G4Indoorv1_argparse.h"
#include "bs_tracing.h"
#include "bs_oswrap.h"
#include "bs_cmd_line.h"

static char library_name[] = "2G4 Indoors channel model";

void component_print_post_help() {
  fprintf(stdout,
    "For more information on the format of the <distance_matrix_file> and the\n"
    "<distance_file> please check: doc/README.txt\n"
    "Note that in utilities/ you can find a MATLAB script to generate a\n"
    "distance file from devices coordinates\n"
    "\n"
    "Note that parameters are processed from left to right and can be overriden\n"
    "\n"
    "Model of the 2.4GHz indoors channel. It models the effects of the multipath\n"
    "fading of the channel.\n"
    "It creates NxN independent channels/statistical processes, all of them with\n"
    "the same basic parameters as specified in the command line\n"
  );
}

ch_2G4I_args_t *args_g;

void cmd_att_found(char * argv, int offset) {
  if ((args_g->attenuation < -100) || (args_g->attenuation > 100)) {
    bs_trace_error("channel: cmdarg: attenuation can only be between -100 and 100dB (%lf)\n",
                   args_g->attenuation);
  }
}
void cmd_attextra_found(char * argv, int offset) {
  if ((args_g->atxtra < -100) || (args_g->atxtra > 100)) {
    bs_trace_error("channel: cmdarg: extra attenuation can only be between -100 and 100dB (%lf)\n",
                   args_g->atxtra);
  }
}
void cmd_exp_found(char * argv, int offset) {
  if ((args_g->distance_exp < 1) || (args_g->distance_exp > 4)) {
    bs_trace_error("channel: cmdarg: distance_exp can only be between 1 and 4 (%lf) '%s'\n",
                   args_g->distance_exp, argv);
  }
}
void cmd_ds_found(char * argv, int offset) {
  if ((args_g->RMS_DelaySpread < 0) || (args_g->RMS_DelaySpread > 100)) {
    bs_trace_error("channel: cmdarg: delay spread can only be between 0 and 100ns (really 0..30ns is the expected range) (%lf) '%s'\n",
                   args_g->RMS_DelaySpread, argv);
  }
  args_g->RMS_DelaySpread = args_g->RMS_DelaySpread*1e-9;
}
void cmd_speed_found(char * argv, int offset) {
  if ((args_g->DopplerSpeed < 0) || (args_g->DopplerSpeed > 10)) {
    bs_trace_error("channel: cmdarg: Doppler speed can only be between 0 and 10m/s (really 0..1.1m/s is the expected range) (%lf) '%s'\n",
                   args_g->DopplerSpeed, argv);
  }
}
void cmd_RiceK_found(char * argv, int offset) {
  if ((args_g->Rice_K < 0) || (args_g->Rice_K > 20)) {
    bs_trace_error("channel: cmdarg: Rice K factor can only be between 0 and 20 (really 0..4 is the expected range) (%lf) '%s'\n",
                   args_g->Rice_K, argv);
  }
}

static char *preset;
void cmd_Preset_found(char * argv, int offset) {
  if (strcasecmp(preset,"Small2") == 0) {
    args_g->RMS_DelaySpread = 10e-9;
    args_g->Rice_K = 0.65;
  } else if (strcasecmp(preset,"Big4") == 0) {
    args_g->RMS_DelaySpread = 15e-9;
    args_g->Rice_K = 0.75;
  } else if (strcasecmp(preset,"Huge3") == 0) {
    args_g->RMS_DelaySpread = 25e-9;
    args_g->Rice_K = 1.0;
  } else if (strcasecmp(preset,"Huge10") == 0) {
    args_g->RMS_DelaySpread = 25e-9;
    args_g->Rice_K = 0.0;
  } else {
    bs_trace_error_line("channel: Preset \"%s\" is not known\n",preset);
  }
}
/**
 * Check the arguments provided in the command line: set args based on it
 * or defaults, and check they are correct
 */
void channel_2G4I_argparse(int argc, char *argv[], ch_2G4I_args_t *args)
{
  args_g = args;
  bs_args_struct_t args_struct[] = {
    /*manual,mandatory,switch,option,   name ,      type,   destination,               callback,      , description*/
    { false, false, false, "at"    ,"att",           'f', (void*)&args->attenuation, cmd_att_found, "default path loss attenuation in dB, used in all NxN paths if <dist_matrix_file> is not provided (default 51.8dB = 4meters) (can be be < 0 for debugging purposes)"},
    { false, false, false, "atextra","atextra",      'f', (void*)&args->atxtra,      cmd_attextra_found, "Extra attenuation to be added to all paths (even if <dist_matrix_file> is provided)"},
    { false, false, false, "dist","dist_matrix_file",'s', (void*)&args->matrix_file_name, NULL, "File containing the distance for each NxN path"},
    { false, false, false, "preset","PreSetName",'s', (void*)&preset, cmd_Preset_found, "One of the predefined presets can be chosen: "
       "{Small2, (Big4), Huge3, Huge10} "
       "This will set the <delay_spread> and <RiceK> as: "
       "Small room (~10m2) where devices are ~2m apart (delay_spread=10ns, riceK=0.65); "
       "Big room (~35m2) where devices are ~4m apart (delay_spread=15ns, riceK=0.75); "
       "Huge room (~100m2) where devices are ~3m apart (delay_spread=25ns, riceK=1); "
       "Huge room (~100m2) where devices are ~10m apart (delay_spread=25ns, riceK=0)"},
    { false, false, false, "speed","doppler_speed",'f', (void*)&args->DopplerSpeed, cmd_speed_found, "Speed of the transmitter, receiver and/or enviroment in m/s (default 1.1m/s = 4km/h = slow walk)"},
    { false, false, false, "ds",  "delay_spread",'f', (void*)&args->RMS_DelaySpread, cmd_ds_found, "RMS delay spread (S) of the channel (in ns) (default 15ns)"},
    { false, false, false, "RiceK", "riceK",'f', (void*)&args->Rice_K, cmd_RiceK_found, "Ricean K factor: ratio of the line of sight power to the scattered power. (Default 2), expected range from 0 to 3 = No line of sight to approx 1.5m distance between tx and rx in a small room"},
    { false, false, false, "exp", "distance_exp",'f', (void*)&args->distance_exp, cmd_exp_found, "Distance exponent for the path loss (default 2)"},
    ARG_TABLE_ENDMARKER
  };

  //set defaults: (walking slowly in a medium room, with devices at around 4 meters from each other)
  args->DopplerSpeed = 1.1;
  args->RMS_DelaySpread = 15e-9;
  args->Rice_K = 0.75;
  args->attenuation = 51.8; //4meters attenuation
  args->atxtra = 0;
  args->matrix_file_name = NULL;
  args->distance_exp = 2.0;

  char trace_prefix[50]; //it will not be used as soon as we get out of this function
  snprintf(trace_prefix,50, "channel: (Indoor) ");
  bs_args_override_exe_name(library_name);
  bs_args_set_trace_prefix(trace_prefix);
  bs_args_parse_all_cmd_line(argc, argv, args_struct);
}
