#!/bin/bash
# Copyright 2018 Oticon A/S
# SPDX-License-Identifier: Apache-2.0
# This test requires the results to be inspected manually

SIMULATION_ID="silly"
VERBOSITY_LEVEL=4
PROCESS_IDS=""; EXIT_CODE=0

function Execute(){
  if [ ! -f $1 ]; then
    echo -e "  \e[91m`pwd`/`basename $1` cannot be found (did you forget to\
 compile it?)\e[39m"
    exit 1
  fi
  timeout 5 $@ & PROCESS_IDS="$PROCESS_IDS $!"
}

BSIM_OUT_PATH="${BSIM_OUT_PATH:-../../../}"
cd ${BSIM_OUT_PATH}/bin

Execute ./bs_device_2G4_playback \
  -v=${VERBOSITY_LEVEL} -s=${SIMULATION_ID} -d=0 -inputf=../components/ext_2G4_channel_Indoorv1/test/data2/0

Execute ./bs_device_2G4_playback \
  -v=${VERBOSITY_LEVEL} -s=${SIMULATION_ID} -d=1 -inputf=../components/ext_2G4_channel_Indoorv1/test/data2/1

Execute ./bs_device_2G4_playback \
  -v=${VERBOSITY_LEVEL} -s=${SIMULATION_ID} -d=2 -inputf=../components/ext_2G4_channel_Indoorv1/test/data2/2

Execute ./bs_2G4_phy_v1 -v=${VERBOSITY_LEVEL} -s=${SIMULATION_ID} -channel=Indoorv1 \
  -argschannel -at=30 -dist=../components/ext_2G4_channel_Indoorv1/test/data2/silly.matrix -atextra=10  -argsmain \
  -D=3 -sim_length=20e6 $@

for PROCESS_ID in $PROCESS_IDS; do
  wait $PROCESS_ID || let "EXIT_CODE=$?"
done
exit $EXIT_CODE #the last exit code != 0
