/*
 * Copyright 2018 Oticon A/S
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef BS_CHANNEL_2G4INDOORS_ARGS_H
#define BS_CHANNEL_2G4INDOORS_ARGS_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct{
  double attenuation;
  double atxtra;
  char *matrix_file_name;
  double distance_exp;
  double RMS_DelaySpread;
  double DopplerSpeed;
  double Rice_K;
} ch_2G4I_args_t;

void channel_2G4I_argparse(int argc, char *argv[], ch_2G4I_args_t *args);

#ifdef __cplusplus
}
#endif

#endif
