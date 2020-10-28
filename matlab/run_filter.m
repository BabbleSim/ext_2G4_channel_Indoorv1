function [output, FilterMemory] = run_filter(input, filtercoefs, FilterMemory)
%
% Copyright 2014 Oticon A/S
% SPDX-License-Identifier: Apache-2.0

FilterMemory = circshift(FilterMemory',1)';
FilterMemory(:,1) = input;

output = FilterMemory*filtercoefs';
