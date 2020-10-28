% Copyright 2014 Oticon A/S
% SPDX-License-Identifier: Apache-2.0


opengl software
set(gcf, 'Renderer', 'opengl')


%%

SimTime = 1; %in seconds
paint = 1;

RMS_DelaySpread = 30e-9; %ns

DopplerSpeed = 4/3.6; %4km/h

RNGSeed = 1341234;
rng(RNGSeed)

fs = 80e6; %channel BW (just the ISM BW)
Nsamples_ChannelAnalysis = 160; %FFT length for the 80MHz
Oversampling_recalc = 16; %how often in time do we recalcultate new taps compared to the Nyquist limit (for that given coherence time/DopplerSpread)
Rice_K = 0; %For LOS: 0..2 tipical for LOS indoors "small" rooms
              %for NLOS a tipical Rice_K = 0, and shadowing therefore of approx. 3.8dB

[ Taps , FFChannelResp, Sim_times ] = Calculate_Fast_Fading(RMS_DelaySpread, DopplerSpeed, Rice_K,  SimTime, Nsamples_ChannelAnalysis, Oversampling_recalc);
FFChannelRespdB = 20*log10(abs(FFChannelResp));

if paint,
    figure1 = figure();
    axes1 = axes('Parent',figure1);
    [x,y] = meshgrid( (0:1/Nsamples_ChannelAnalysis:1-1/Nsamples_ChannelAnalysis)*fs/1e6 , Sim_times);
    mesh(x,y,FFChannelRespdB ,'Parent',axes1)
    xlabel('Freq MHz','Parent',axes1);
    ylabel('time seconds','Parent',axes1);
    zlabel('Fading dB','Parent',axes1);
    view(axes1,[39.5 48]);
    title({['ISM band fast fading. DopplerSpeed = ' num2str(round(DopplerSpeed*10)/10) 'm/s'] , ['RMS DelaySpread = ' num2str(RMS_DelaySpread*1e9) 'ns ; K = ' num2str(Rice_K)]},'Parent',axes1);
end

%figure(); plot(Sim_times,  FFChannelResp(:,1)); xlabel('seconds'); ylabel('dBm');
%figure(); plot( (0:1/Nsamples_ChannelAnalysis:1-1/Nsamples_ChannelAnalysis)*fs/1e6, FFChannelResp(1,:)); xlabel('MHz'); ylabel('dBm');


%% Calculate (roughly) the equivalent loss the receiver sees and the ISI due to fading
%we assume an "ideal" passband filter of the BW below

BW = 3; %MHz
ChannelCenters = 2:2:80-BW/2;

FTapSize = 80/Nsamples_ChannelAnalysis;

ISI_SNR = zeros(length(Sim_times),length(ChannelCenters));
AverageFadeLevel = zeros(length(Sim_times),length(ChannelCenters));

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
  AverageFadeLevel(:, ChannelCenter == ChannelCenters) = 10*log10(D);
end

if paint,
    figure1 = figure();
    axes1 = axes('Parent',figure1);
    [x,y] = meshgrid( ChannelCenters , Sim_times);
    contourf(x,y,max(AverageFadeLevel,-35),'Parent',axes1)
    xlabel('Channel','Parent',axes1);
    ylabel('time seconds','Parent',axes1);
    zlabel('Effective fading dB','Parent',axes1);
    %view(axes1,[39.5 48]);
    c = colorbar
    ylabel(c,'Effective fading dB');
    
    title({['Effective fading, BW ' num2str(BW) 'MHz'];['ISM band fast fading. DopplerSpeed = ' num2str(round(DopplerSpeed*10)/10) 'm/s ; RMS DelaySpread = ' num2str(RMS_DelaySpread*1e9) 'ns ; Rice K = ' num2str(Rice_K)]},'Parent',axes1);

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

return

%% Trial to understand the ISI we face at different BWs in each channel
%cell just for understanding/investigating purposes (expanding/testing the calculation done in the previous cell)

Channel = 0; %0..39
BW = 1; %MHz

TimeChops = 1:length(Sim_times);
for TimeChop = TimeChops,
    Fcenter = Channel*2 + 2;
    Fstart = Fcenter - BW/2;
    Fend = Fcenter + BW/2;
    FTapSize = 80/Nsamples_ChannelAnalysis;
    FirstFTap = Fstart/FTapSize;
    LastFTap = Fend/FTapSize;

    ChannelResp = FFChannelResp( TimeChop, FirstFTap+1:LastFTap ); 
    
    %%thoughts
    N= length(ChannelResp);
    D = sum(ChannelResp)/N;
    Dd = 20*log10(abs(D));
    AllOther = sqrt(sum(abs(ChannelResp).^2)/N - abs(D)^2);
    %AllOtherdB = 20*log10(AllOther);
    %fprintf('Freq: DcPower = %2.3f  ; AllOtherPower = %.3f   \n', abs(D)^2, AllOther^2);
    %fprintf('Freq: DcPowdB = %2.3fdB; AllOtherPowdB = %.3fdB \n',  Dd, AllOtherdB);
    %%endofthoughts
    
    figure(1);
    plot(Fstart+FTapSize:FTapSize:Fend , 20*log10(abs(ChannelResp./max(ChannelResp))));
    xlabel('MHz + 2400');
    ylabel('dB');
    title([num2str(20*log10(max(abs(ChannelResp)))) 'dBs']);
        
    %equivalent to cutting with a rectangular window in freq == like if the Rx filter was a perfect sync
    TimeResp = fftshift(ifft(ChannelResp));
    
    %%thoughts
    AllTime = sum(abs(TimeResp).^2);
    AllOtherTime = sum(abs(TimeResp).^2) - max(abs(TimeResp))^2;
    MaxTap= max(abs(TimeResp))^2;
    values = sort(20*log10(abs(TimeResp./max(TimeResp))),'descend');
    valuesStore(TimeChop,:) = values(2:11)-values(2);
    %fprintf('Time: MaxTapPower = %2.3f  ; AllOtherPower = %.3f  ; All   = %.3f \n', MaxTap, AllOtherTime, AllTime);
    %fprintf('Time: MaxTapPowdB = %2.3fdB; AllOtherPowdB = %.3fdB; AlldB = %.3fdB \n\n',10*log10(MaxTap), 10*log10(AllOtherTime), 10*log10(AllTime));
    %%end of thoughts
    
    figure(2);
    plot(0:1/BW:(1/FTapSize-1/BW) , 20*log10(abs(TimeResp./max(TimeResp))),'o');
    xlabel('us');
    ylabel('dB');
    title([num2str(Sim_times(TimeChop)) 's   ;   ' num2str(20*log10(max(abs(TimeResp)))) 'dB']);
    ylim([-50 0]);
    grid on;
    
    %%more thoughts 
    ISI_SNRLoss = 10*log10(MaxTap/AllOtherTime);
    ISI_SNRLossAprox = 20*log10(D/AllOther);
    Fading_LevelLoss = Dd;
    fprintf('Fading_LevelGain = %.3fdB ; ISI = %.3fdB ; ISI_Aprox = %.3fdB \n',Fading_LevelLoss, -ISI_SNRLoss, -ISI_SNRLossAprox);
        
    %%end of more thoughts
    %pause
end


return 


%% Equivalent sync we filter with by just picking the ChannelResp (== ideal rectangular filter in freq)
%cell just for understanding/investigating purposes

BW = 4;
Fstart = 0 - BW/2;
Fend = 0+ BW/2;
FTapSize = 80/Nsamples_ChannelAnalysis;
FirstFTap = Fstart/FTapSize;
LastFTap = Fend/FTapSize;
ChannelFilter = [ zeros(1,FirstFTap) ones(1,length(FirstFTap+1:LastFTap)) zeros(1,Nsamples_ChannelAnalysis-LastFTap) ];
ChannelFilter = [ ones( 1,LastFTap - max(FirstFTap,0) )   zeros(1, Nsamples_ChannelAnalysis - (LastFTap-FirstFTap) )  ones(1,max(-FirstFTap,0)) ];
TimeResp = fftshift(ifft(ChannelFilter));
figure(2);
plot((0:1/80:(1/FTapSize-1/80)) , (abs(TimeResp./max(TimeResp))));
xlabel('us');
ylabel('n.u.');
grid on;

return

