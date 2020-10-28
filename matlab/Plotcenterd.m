function [] = Plotcenterd(signal, Ts, varargin)
%
% Copyright 2014 Oticon A/S
% SPDX-License-Identifier: Apache-2.0

center = find(signal == max(signal));

x_axis = Ts*( (1:length(signal)) - center);
plot(x_axis, 20*log10(abs(signal./max(signal))),varargin{:});
