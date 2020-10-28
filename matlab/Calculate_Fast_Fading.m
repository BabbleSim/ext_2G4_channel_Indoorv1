function [Taps, ChannelResp, Sim_times] = Calculate_Fast_Fading(RMS_DelaySpread, DopplerSpeed, Rice_K, SimTime, Nsamples_ChannelAnalysis, Oversampling_recalc)
% function [Taps, ChannelResp, Sim_times] = Calculate_Fast_Fading(RMS_DelaySpread, DopplerSpeed, Rice_K, SimTime, Nsamples_ChannelAnalysis, Oversampling_recalc)
% 
% Calculate a fast fading process for the ISM band with a given RMS delay spread
% and doppler spread (speed relative to surroundings or of objects movign around)
% it assumes a flat doppler spectrum, and a exponential power delay profile
%
% Note that this fast fading process has NO MEAN GAIN
% 
% inputs
%  RMS_DelaySpread in seconds (e.g.: 30e-9)
%     rms delay corresponds to a e^(-1) decay in the power delay
%     profile == -8.68dB. A 30dB decay == takes 3.45 times longer
%     86 ns == ITU recommended profile for WLAN indoors.
%     indoor residential has a delay spread in the range of 10 to 30ns
%     big open areas up to 170ns in extreme cases (see ITU P 1238, table 5)
%
%  DopplerSpeed in m/s (speed of movement relative to the transmitter of
%  the receiver or the surroundings)
%
%  Rice_K: Rice K factor: the ratio between the power in the direct path
%             and the power in the other, scattered, paths K = v^2 / (2*sigma^2)
%             for 2.4GHz, indoors, and "small" rooms
%             tipical values for LOS = 0..2 tipical for LOS indoors "small" rooms (extremes of maybe 3 or 4)
%                  (the understanding is that as you get far from the transmitter the direct path fades faster than the scatter)
%                   see LustmannPorrat_2010_Kfactor.pdf
%             tipical value for NLOS = 0 (remember to add then an equivalent shadow
%             factor to the output of 3.8dB if NLOS)
%             The total Rice power is
%             1 = omega =(1+K)*2*sigma^2 where sigma is the std of each component
%             G(0,sigma) ( http://en.wikipedia.org/wiki/Rician_fading )
%
%  SimTime in seconds
% 
%  Nsamples_ChannelAnalysis : How many samples for the channel response
%  (divisions of the 80MHz band). Recommended 128
% 
%  Oversampling_recalc: how often in time do we recalcultate new taps
%     compared to the Nyquist limit (for that given coherence time/DopplerSpread)
%     recommended _16_
%
%
% For a justification of the exponential power delay profile assumption see for ex:
% multipath channel parameters for the indoor ratio at 2.4GHz ISM band.pdf
% the ITU recommendation P.1238-7
% ITU model for indoor attenuation_R-REC-P.1238-7-201202-I!!PDF-E.pdf
% or the proposed channel for indoors in 
% 15-12-0459-07-0008-tg8-channel-models.doc
% (exponential)
%
% Copyright 2014 Oticon A/S
% SPDX-License-Identifier: Apache-2.0


MaxDopplerSpread = DopplerSpeed/299.792458e6*2.4e9; %(approx 9Hz at 4km/h)
Ts = 1/160e6; %6.25ns: enough to get the 80MHz BW (160MHz sampling)
Nsamples_ChannelAnalysis = Nsamples_ChannelAnalysis*2; %FFT length for the _160MHz_

%we normalize the total Rice power to 1 = omega
taps_gain = sqrt( 1 / ( 2* ( 1 + Rice_K ) ) ); %== gausians std (not std^2)
rice_offset = sqrt(Rice_K / (1+Rice_K)); %in amplitude

[taps_std, ~] = Taps_level_generator(RMS_DelaySpread, Ts); %generate Rayleight taps std level and offsets
n_taps = length(taps_std);

if MaxDopplerSpread ~= 0,
    Recalc_time = 1/MaxDopplerSpread/Oversampling_recalc;
else
    Recalc_time = 1e6;
end

[filter_coefs, FilterCorrection] = Calculate_doppler_filter(Oversampling_recalc); %Generate low pass filter for the Rayleight generators (flat doppler spread)
RandomGain = taps_gain*FilterCorrection;

Init = length(filter_coefs)-1;
FilterMemory = zeros(n_taps,length(filter_coefs));

Taps          = zeros(floor(SimTime/Recalc_time),n_taps);
ChannelResp   = zeros(floor(SimTime/Recalc_time),Nsamples_ChannelAnalysis/2);
Sim_times     = zeros(1,floor(SimTime/Recalc_time));
index = 0;
for time = -Init*Recalc_time:Recalc_time:max(SimTime, Recalc_time),
    Taps_Rnd = ( randn(1,n_taps) + 1j*randn(1,n_taps))*RandomGain; %Generate random input to Rayleight generators
    [Tap_filtered, FilterMemory]  = run_filter(Taps_Rnd,filter_coefs,FilterMemory); %filter to account for Doppler spread
    %each Tap_filtered has std and std = 1 at this point
    if time < 0,%the first Init rounds are just to ramp up the filters
        continue;
    end
    index = index + 1;
    Tap_filtered = Tap_filtered.*taps_std'; %we scale them with their relative level
    Tap_filtered(1) = Tap_filtered(1) + rice_offset;
    Taps(index,:) = Tap_filtered; %we store them for later debugging
    
    ChannelResponse = fft([Tap_filtered ; zeros(Nsamples_ChannelAnalysis-length(Tap_filtered),1) ] ); %FFT resolution of 160MHz/Nsamples_ChannelAnalysis and span of +-80MHz (fs =160MHz)
    ChannelResp(index,:) = ChannelResponse(1:Nsamples_ChannelAnalysis/2); %we keep the first 80MHz
    Sim_times(index) = time;
end

