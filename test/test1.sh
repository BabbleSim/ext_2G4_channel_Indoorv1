#!/usr/bin/env bash
# Copyright 2018 Oticon A/S
# SPDX-License-Identifier: Apache-2.0

SIMULATION_ID="Silly"
VERBOSITY_LEVEL=5
PROCESS_IDS=""; EXIT_CODE=0
function Execute(){
  timeout 5 $@ & PROCESS_IDS="$PROCESS_IDS $!"
}

BSIM_OUT_PATH="${BSIM_OUT_PATH:-../../../}"
cd ${BSIM_OUT_PATH}/bin

Execute ./bs_2G4_phy_v1 -v=${VERBOSITY_LEVEL} -s=${SIMULATION_ID} -D=2 -sim_length=5e5 $@ \
-channel=Indoorv1 -argschannel -RiceK=1.5 -speed=0.5 -preset=Huge3 -at=-55 -ds=30 \
-exp=2.0 -dist=../components/ext_2G4_channel_Indoorv1/test/data1/trial_case_1.matrix 

Execute ./bs_device_2G4_playback -d=0 -v=${VERBOSITY_LEVEL} -s=${SIMULATION_ID} \
-inputf=../components/ext_2G4_channel_Indoorv1/test/data1/Play1.d_0 

Execute ./bs_device_2G4_playback -d=1 -v=${VERBOSITY_LEVEL} -s=${SIMULATION_ID} \
-inputf=../components/ext_2G4_channel_Indoorv1/test/data1/Play1.d_1 

for PROCESS_ID in $PROCESS_IDS; do
  wait $PROCESS_ID || let "EXIT_CODE=$?" 
done
exit $EXIT_CODE #the last exit code != 0
