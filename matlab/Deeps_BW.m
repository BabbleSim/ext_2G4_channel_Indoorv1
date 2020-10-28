%Calculate the CDF of the BW of the dips of the fast fading for a given
%depth
%
% Copyright 2014 Oticon A/S
% SPDX-License-Identifier: Apache-2.0



SimTime = 2000; %in seconds

RMS_DelaySpread = 70e-9; %indoor residential with delay spread in the range of 20 to 150ns
                  %big open areas up to 500ns in extreme cases
            %rms delay corresponds to a e^(-1) decay in the power delay
            %profile == -8.68dB. A 30dB decay == takes 3.45 times longer
            %86 ns == ITU default

DopplerSpeed = 4/3.6; %4km/h

RNGSeed = 1341234;
rng(RNGSeed)

fs = 80e6; %channel BW (just the ISM BW)
Nsamples_ChannelAnalysis = 512; %FFT length for the 80MHz
Oversampling_recalc = 3; %how often in time do we recalcultate new taps compared to the Nyquist limit (for that given coherence time/DopplerSpread)
Rice_K = 0; %For LOS: 0..2 tipical for LOS indoors "small" rooms   
              %for NLOS a tipical Rice_K = 0, 

[ Taps , FFChannelResp, Sim_times ] = Calculate_Fast_Fading(RMS_DelaySpread, DopplerSpeed, Rice_K,  SimTime, Nsamples_ChannelAnalysis, Oversampling_recalc);
FFChannelRespdB = 20*log10(abs(FFChannelResp));



%% calculate dips CDFs
depths = [-10 -15 -20 -25]; %dBs

tap_thicknes = fs/Nsamples_ChannelAnalysis;
range = 0:1:ceil(12e6/(fs/Nsamples_ChannelAnalysis))-1 ; %the range of the data, with enough granularity
CDF = zeros(length(depths), length(range) );
pdfs = zeros(length(depths), length(range) );
for Depth_Th = depths,
    UnderTh = FFChannelRespdB < Depth_Th;
    found = 0;
    Durations = [];
    for i = 1:size(UnderTh,1),
        %UnderTh(i,:); %cut over BW in a given moment
        dsig = diff([ 0 UnderTh(i,:) 0 ]);
        startIndexes = find(dsig > 0);
        endIndexes = find(dsig < 0)-1;
        duration = endIndexes-startIndexes+1;
        Durations(end+1 : (end+length(duration)) ) = duration;
    end

    pdf = hist(Durations,range);
    CDF(Depth_Th == depths,:) = cumsum(pdf)/length(Durations);
    pdfs(Depth_Th == depths,:) = pdf/length(Durations);
end

figure(10); hold off;
plot(range*tap_thicknes/1e6,CDF);
title(['CDF of the BW of a fading dip of X dBs for a r.m.s. delay spread of ' num2str(round(RMS_DelaySpread*1e9*100)/100) 'ns']);
ylabel('prob');
xlabel('MHz');
grid on
Names = [num2str(depths') ,  repmat('dBs',length(depths),1)];
legend(Names,'Location','Best');

figure(11); hold off;
plot(range*tap_thicknes/1e6,pdfs);
title(['pdf of the BW of a fading dip of X dBs for a r.m.s. delay spread of ' num2str(round(RMS_DelaySpread*1e9*100)/100) 'ns']);
ylabel('prob');
xlabel('MHz');
grid on
Names = [num2str(depths') ,  repmat('dBs',length(depths),1)];
legend(Names,'Location','Best');

