function [taps, times] = Taps_level_generator(RMS_DelaySpread, Ts)
%function [taps, times] = Taps_level_generator(RMS_DelaySpread, Ts)
%
% Generate taps levels corresponding to a exponential delay spread
% with rms delay spread = RMS_DelaySpread 
% and sampled at Ts
%
% taps are rayleigh taps levels in natural units (amplitude = sigma)
% times are the times of those taps
% they are normalized to have power 1 all together
%
% Copyright 2014 Oticon A/S
% SPDX-License-Identifier: Apache-2.0


MaxT = 3.5*RMS_DelaySpread; %downto -30dB
times = 0:Ts:MaxT;

taps = exp(-times/RMS_DelaySpread);

taps = taps/sqrt(sum(taps.^2));
