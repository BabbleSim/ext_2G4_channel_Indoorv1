# Copyright 2018 Oticon A/S
# SPDX-License-Identifier: Apache-2.0

%this file specifies the devices  positions over time into the devices{} cell array
% this is used by create distance_files.m to generate the 2G4_channel_Indoorv1 distance files

%sampling of the position of the device x over time: (time shall be in ascending order)
device{1}.x = [1 1 10];
device{1}.y = [0 0 0];
device{1}.time = [0 30 100];

device{2}.x = 3;
device{2}.y = 0;
device{2}.time = 0;

device{3}.x = [2 2];
device{3}.y = [2 12]; %1 m/s going out in y axes
device{3}.time = [0 100];

device{4} = device{3};
device{4}.x = device{4}.x + 1; 

outputfileprefix = 'trial_case_2';