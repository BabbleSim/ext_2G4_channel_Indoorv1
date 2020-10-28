% Copyright 2014 Oticon A/S
% SPDX-License-Identifier: Apache-2.0


opengl software
set(gcf, 'Renderer', 'opengl')


%%

SimTime = 1; %in seconds
paint = 1;

PowerTx = 20; %dBm
ReceiverGain = -15; %dB

RMS_DelaySpread = 30e-9; 

DopplerSpeed = 4/3.6; %4km/h

RNGSeed = 1341234;
rng(RNGSeed)

fs = 80e6; %channel BW (just the ISM BW)
Nsamples_ChannelAnalysis = 256; %FFT length for the 80MHz
Oversampling_recalc = 16; %how often in time do we recalcultate new taps compared to the Nyquist limit (for that given coherence time/DopplerSpread)
Rice_K = 0; %For LOS: 0..2 tipical for LOS indoors "small" rooms
              %for NLOS a tipical Rice_K = 0, and shadowing therefore of approx. 3.8dB

[ Taps , FFChannelResp, Sim_times ] = Calculate_Fast_Fading(RMS_DelaySpread, DopplerSpeed, Rice_K,  SimTime, Nsamples_ChannelAnalysis, Oversampling_recalc);
FFChannelRespdB = 20*log10(abs(FFChannelResp));

WalkSpeed = 0/3.6; %4km/h
distance = 5 + (Sim_times*WalkSpeed); %meters

PL = PathLoss(distance) + PowerTx + ReceiverGain;
PLv = repmat(PL',1,size(FFChannelRespdB,2));

if paint,
    figure1 = figure();
    axes1 = axes('Parent',figure1);
    [x,y] = meshgrid( (0:1/Nsamples_ChannelAnalysis:1-1/Nsamples_ChannelAnalysis)*fs/1e6 , Sim_times);
    %mesh(x,y,FFChannelRespdB - PLv,'Parent',axes1)
    mesh(x,y,FFChannelRespdB ,'Parent',axes1)
    %surf(x,y,FFChannelRespdB ,'Parent',axes1,'LineStyle','none'); colorbar('peer',axes1);
    xlabel('Freq MHz','Parent',axes1);
    ylabel('time seconds','Parent',axes1);
    zlabel('Fading dB','Parent',axes1);
    view(axes1,[39.5 48]);
    title(['ISM band fast fading. DopplerSpeed = ' num2str(round(DopplerSpeed*10)/10) 'm/s ; RMS DelaySpread = ' num2str(RMS_DelaySpread*1e9) 'ns ; K = ' num2str(Rice_K)],'Parent',axes1);
end

%figure(); plot(Sim_times,  FFChannelResp(:,1)); xlabel('seconds'); ylabel('dBm');
%figure(); plot( (0:1/Nsamples_ChannelAnalysis:1-1/Nsamples_ChannelAnalysis)*fs/1e6, FFChannelRespdB(1,:)); xlabel('MHz'); ylabel('dBm');
