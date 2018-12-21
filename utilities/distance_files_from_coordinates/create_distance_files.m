# Copyright 2018 Oticon A/S
# SPDX-License-Identifier: Apache-2.0

%Script to generate the distance files for the 2G4_channel_Indoorv1 from a
%matlab file which specified the devices coordinates over time

%TOOPT: could make sense to rewrite this in Python to save the matlab license

clear
OutputDir = ['..' filesep '..' filesep 'test' filesep 'data' filesep]; %directory relative to where this script is executed (where the result of this script is stored)
RelativePath = '../components/2G4_channel_Indoorv1/test/data1/'; %path from bin/ into the place where the files are stored (this path is to be used while reading)

%we load the case (number of devices as positions over time):
trial_case_1  %<---------- here is where the define the distances over time

generate_distance_files(device, outputfileprefix, OutputDir, RelativePath);
