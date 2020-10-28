function [TimeResp, Times, FreqResp, Freqs] = TimeResp_Gausian(SymbolFs, BT, Ts, Duration, Fs, MaxFreq)
%function [TimeResp, Times, FreqResp, Freqs] = TimeResp_Gausian(SymbolFs, BT, Ts, Duration, Fs, MaxFreq)
% generate the impulse response (and freq. resp.) of a gausian pulse shaping filter
% and the freq. response
%
% SymbolFs: data rate (e.g. = 1e6)
% BT : Gaussian BT (B*SymbolTs, where B = 3dB BW of the filter (between DC and B), SymbolTs = 1/SymbolFs) (e.g. 0.5)
%
% For the time response:
% Duration: impulse response duration (e.g. 8*1/SymbolFs)
% Ts: impulse response sampling resolution (e.g. 1/80 * 1/SymbolFs)
%
% For the frequency response: Note that this is the response of a *continous* exponential, the produced exponential will have a different response as it is not infinitely oversampled == its response is aliased == flatter.
% Fs: frequency response sampling resolution (e.g. SymbolFs/10)
% MaxFreq: the response will be generated in -MaxFreq:Fs:MaxFreq
%
% Copyright 2014 Oticon A/S
% SPDX-License-Identifier: Apache-2.0

%http://www.mathworks.se/help/signal/examples/fir-gaussian-pulse-shaping-filter-design.html

SymbolTs = 1/SymbolFs; %Ts of the actual data symbol

B = BT/SymbolTs; %3dB BW of the filter = �B (between DC and B)
a= 1/BT*sqrt(log(2)/2);

%Hideal = exp(-a^2*(f/fs).^2);
%h = sqrt(pi)/a * exp(- (pi*t/a/Ts).^2)
%where a = 1/(B*Ts) * sqrt(ln(2)/2)

Times = -Duration/2:Ts:Duration/2-Ts;
TimeResp = exp(-(pi*Times/a/SymbolTs).^2); % factor sqrt(pi)/a removed to have max amplitude 1
Freqs = -MaxFreq:Fs:MaxFreq;
FreqResp = exp(-(a*(Freqs/SymbolFs)).^2);
