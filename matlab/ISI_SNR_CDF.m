function [CDF, bins, pdf] = ISI_SNR_CDF(BW, Rice_K, RMS_DelaySpread, DopplerSpeed, RNGSeed, SimTime, dBRange)
%function [CDF, bins, pdf] = ISI_SNR_CDF(BW, Rice_K, RMS_DelaySpread, DopplerSpeed, RNGSeed, SimTime, dBRange);
% BW = 3; %MHz
% SimTime = 1; %in seconds
% RMS_DelaySpread = 10e-9; %ns
% DopplerSpeed = 4/3.6; %4km/h
% RNGSeed = 1341234;
% Rice_K = 0; %For LOS: 0..2 tipical for LOS indoors "small" rooms
%             %for NLOS a tipical Rice_K = 0, and shadowing therefore of approx. 3.8dB
%
% Copyright 2014 Oticon A/S
% SPDX-License-Identifier: Apache-2.0



paint = 0;
paint_channel = 0;
paint_countours = 0;

%%

rng(RNGSeed)

fs = 80e6; %channel BW (just the ISM BW)
Nsamples_ChannelAnalysis = 160; %FFT length for the 80MHz
Oversampling_recalc = 16; %how often in time do we recalcultate new taps compared to the Nyquist limit (for that given coherence time/DopplerSpread)
[ ~ , FFChannelResp, Sim_times ] = Calculate_Fast_Fading(RMS_DelaySpread, DopplerSpeed, Rice_K,  SimTime, Nsamples_ChannelAnalysis, Oversampling_recalc);


%% Calculate (roughly) the equivalent loss the receiver sees and the ISI due to fading
%we assume an "ideal" passband filter of the BW below


ChannelCenters = 2:2:80-BW/2;
FTapSize = 80/Nsamples_ChannelAnalysis;

ISI_SNR = zeros(length(Sim_times),length(ChannelCenters));

for ChannelCenter = ChannelCenters,
  Fstarts = ChannelCenter - BW/2;
  Fends = ChannelCenter + BW/2;
  FirstFTap = Fstarts/FTapSize;
  LastFTap = Fends/FTapSize;
  ChannelResp = FFChannelResp( :, FirstFTap+1:LastFTap);
  N= LastFTap - FirstFTap;
  D = abs(sum(ChannelResp,2)/N).^2; %tap 0 power
  AllOther = sum(abs(ChannelResp).^2,2)/N - D; %all others taps powers
  ISI_SNR(:, ChannelCenter == ChannelCenters) = 10*log10(D./AllOther);
end

if paint_countours,
    figure1 = figure();
    axes1 = axes('Parent',figure1);
    [x,y] = meshgrid( ChannelCenters , Sim_times);
    %mesh(x,y,min(ISI,0) ,'Parent',axes1)
    contourf(x,y,min(max(ISI_SNR,0),34) ,'Parent',axes1)
    xlabel('Channel','Parent',axes1);
    ylabel('time seconds','Parent',axes1);
    c= colorbar
    ylabel(c,'SNR dB');
    %view(axes1,[39.5 48]);
    title({['ISI SNR due to Multipath, BW ' num2str(BW) 'MHz'];['ISM band fast fading. DopplerSpeed = ' num2str(round(DopplerSpeed*10)/10) 'm/s ; RMS DelaySpread = ' num2str(RMS_DelaySpread*1e9) 'ns ; Rice K = ' num2str(Rice_K)]},'Parent',axes1);
end

%% Calculate probability of ISI SNR
[pdf, bins] = hist(reshape(ISI_SNR,1,[]),dBRange);
%hold off; plot(bins,pdf);
CDF = (cumsum(pdf)/numel(ISI_SNR));

if paint == 1,
  figure(); hold off;
  hplot = plot(bins,CDF,'DisplayName','calc');
  %hold on;
  %plot(10,CDF(bins==10),'ro','DisplayName','marker');
  makedatatip(hplot,find(bins==10,1));
  xlim([0 30]);
  ylim([1e-4 1]);
  grid on;
  ylabel('Probability');
  xlabel('SNR in dBs');
  title({['ISI, probabiltiy of a SNR <= X dBs, BW ' num2str(BW) 'MHz'];['ISM band fast fading. DopplerSpeed = ' num2str(round(DopplerSpeed*10)/10) 'm/s ; RMS DelaySpread = ' num2str(RMS_DelaySpread*1e9) 'ns ; Rice K = ' num2str(Rice_K)]},'Parent',axes1);
end

