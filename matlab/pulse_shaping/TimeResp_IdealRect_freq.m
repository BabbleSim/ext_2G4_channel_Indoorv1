function [TimeResp, Times] = TimeResp_IdealRect_freq(BW, Ts, Duration)
%function [TimeResp, Times] = TimeResp_IdealRect_freq(BW, Ts)
% generate the impulse response of an ideal rectangular filter in freq of response = 1 in �BW/2, and 0 outside
% Ts = 1/80e6
% BW = 1,2,4e6 MHz
% Duration = 640*Ts = 8e-6
%
% Copyright 2014 Oticon A/S
% SPDX-License-Identifier: Apache-2.0


%Ideal rect. window in freq,

fs = 1/Ts;
Nsamples_ChannelAnalysis = ceil(Duration/Ts); %=> fs/Nsamples_ChannelAnalysis == tap size == max duration of channel response
FTapSize = 1/Duration;

Fstart = 0 - BW/2;
Fend = 0+ BW/2;
FirstFTap = Fstart/FTapSize;
LastFTap = Fend/FTapSize;

ChannelFilter = [ ones( 1,LastFTap - max(FirstFTap,0) )   zeros(1, Nsamples_ChannelAnalysis - (LastFTap-FirstFTap) )  ones(1,max(-FirstFTap,0)) ];

TimeResp = fftshift(ifft(ChannelFilter));
Times = (-Duration/2):Ts:(Duration/2-Ts);
