function [pathloss_indB] = PathLoss(distance)
%
% calculate the path loss due to a given distance for 2.4GHz
% distance in meters
%
% Copyright 2014 Oticon A/S
% SPDX-License-Identifier: Apache-2.0


pathloss_indB = 20*log10(2.4e3) + 20*log10(distance) - 28;