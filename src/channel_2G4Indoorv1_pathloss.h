/*
 * Copyright 2018 Oticon A/S
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef BS_channel_2G4_Pathloss_H
#define BS_channel_2G4_Pathloss_H

#include "bs_types.h"
#include "channel_2G4Indoorv1_argparse.h"

#ifdef __cplusplus
extern "C" {
#endif

void InitPathLoss(ch_2G4I_args_t *args, uint n_devs);
double CalculatePathLoss(uint tx, uint rx, bs_time_t now);
void clearPathLoss();

#ifdef __cplusplus
}
#endif

#endif
