% Copyright 2014 Oticon A/S
% SPDX-License-Identifier: Apache-2.0


%The conclusion from this script is that the pulse shaping will only worsen the ISI
% compared to what was calculated in trial_Generate_Multipath_fromFF.m which only accounted for an "ideal" 
% pulse shape/IF filter being a rectangular filter in frequency

%this script is very rough in both its calculations,
%for the frequency version we assumme that the 

%%


opengl software
set(gcf, 'Renderer', 'opengl')
addpath('./pulse_shaping/')

%%

SimTime = 1; %in seconds
paint = 1;

RMS_DelaySpread = 10e-9; 

DopplerSpeed = 4/3.6; %4km/h

RNGSeed = 1341234;
rng(RNGSeed)

fs = 80e6; %channel BW (just the ISM BW)
Nsamples_ChannelAnalysis = 640; %FFT length for the 80MHz
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
    %surf(x,y,FFChannelRespdB ,'Parent',axes1,'LineStyle','none'); colorbar('peer',axes1);
    xlabel('Freq MHz','Parent',axes1);
    ylabel('time seconds','Parent',axes1);
    zlabel('Fading dB','Parent',axes1);
    view(axes1,[39.5 48]);
    title(['ISM band fast fading. DopplerSpeed = ' num2str(round(DopplerSpeed*10)/10) 'm/s ; RMS DelaySpread = ' num2str(RMS_DelaySpread*1e9) 'ns ; K = ' num2str(Rice_K)],'Parent',axes1);
end

%figure(); plot(Sim_times,  FFChannelResp(:,1)); xlabel('seconds'); ylabel('dBm');
%figure(); plot( (0:1/Nsamples_ChannelAnalysis:1-1/Nsamples_ChannelAnalysis)*fs/1e6, FFChannelResp(1,:)); xlabel('MHz'); ylabel('dBm');


%% Calculate (roughly) the equivalent loss the receiver sees and the ISI due to fading
%we assume an "ideal" passband filter of the BW below + 
% a gausian pulse shape (with NO rectangle)

Filter = 'Gaussian';
%Filter = 'None';

SymbolFs = 1e6;
BT = 0.5;
BW_IF = SymbolFs*8; %the modem oversamples by 8 => the pulse shaping has this sampling freq.
Sampling_Fs = SymbolFs;

Ts = 1/80e6;
Duration = Nsamples_ChannelAnalysis/80e6;
FTapSize = 80e6/Nsamples_ChannelAnalysis;
FTapSize_MHz = 80/Nsamples_ChannelAnalysis;

%Gaus filter alone:
[Resp_PulseShape, ~, HGaus, FreqGaus] = TimeResp_Gausian(SymbolFs, BT, 1/BW_IF, Duration, FTapSize, BW_IF/2);
%PulseShapeInFreq = repmat( HGaus(1:end-1), size(FFChannelResp,1),1 ); %this cant be used as it is an analog filter response not accounting for the aliasing of the response

H3Gaus = fftshift(abs(fft(fftshift(Resp_PulseShape)))); H3Gaus = H3Gaus/max(H3Gaus);
DiffTaps = (BW_IF-Sampling_Fs)/FTapSize;
H4Gaus = H3Gaus(DiffTaps/2+(1:end-DiffTaps));
PulseShapeInFreq = repmat( H4Gaus, size(FFChannelResp,1),1 );

%No filter:
if strcmp(Filter ,'None') == 1,
  PulseShapeInFreq = ones(size(PulseShapeInFreq));
end
%

ChannelCenters = 2:2:78;

ISI = zeros(length(Sim_times),length(ChannelCenters));
AverageFadeLevel = zeros(length(Sim_times),length(ChannelCenters));

for ChannelCenter = ChannelCenters,
  Fstarts = ChannelCenter - Sampling_Fs/1e6/2;
  Fends = ChannelCenter + Sampling_Fs/1e6/2;
  FirstFTap = Fstarts/FTapSize_MHz;
  LastFTap = Fends/FTapSize_MHz;
  ChannelResp = FFChannelResp( :, FirstFTap+1:LastFTap).*PulseShapeInFreq; %still we assume that appart from the pulse shaping we have a perfect filtering
  N= LastFTap - FirstFTap;
  D = abs(sum(ChannelResp,2)/N).^2; %tap 0 power
  AllOther = sum(abs(ChannelResp).^2,2)/N - D; %all others taps powers
  ISI(:, ChannelCenter == ChannelCenters) = -10*log10(D./AllOther);
  AverageFadeLevel(:, ChannelCenter == ChannelCenters) = 10*log10(D);
end

if paint,
%     figure1 = figure();
%     axes1 = axes('Parent',figure1);
%     [x,y] = meshgrid( ChannelCenters , Sim_times);
%     contourf(x,y,AverageFadeLevel ,'Parent',axes1)
%     xlabel('Channel','Parent',axes1);
%     ylabel('time seconds','Parent',axes1);
%     %zlabel('Effective fading dB','Parent',axes1);
%     %view(axes1,[39.5 48]);
%     colorbar
%     title({['Effective fading, BW ' num2str(SymbolFs/1e6) ' (' num2str(BW_IF/1e6) ') MHz ; ' Filter];['ISM band fast fading. DopplerSpeed = ' num2str(round(DopplerSpeed*10)/10) 'm/s ; RMS DelaySpread = ' num2str(RMS_DelaySpread*1e9) 'ns ; Rice K = ' num2str(Rice_K)]},'Parent',axes1);

    figure1 = figure();
    axes1 = axes('Parent',figure1);
    [x,y] = meshgrid( ChannelCenters , Sim_times);
    %mesh(x,y,min(ISI,0) ,'Parent',axes1)
    contourf(x,y,max(min(ISI,0),-30) ,'Parent',axes1)
    xlabel('Channel','Parent',axes1);
    ylabel('time seconds','Parent',axes1);
    %zlabel('dB','Parent',axes1);
    colorbar
    %view(axes1,[39.5 48]);
    title({['ISI due to Multipath, BW ' num2str(SymbolFs/1e6) ' (' num2str(BW_IF/1e6) ') MHz ; ' Filter];['ISM band fast fading. DopplerSpeed = ' num2str(round(DopplerSpeed*10)/10) 'm/s ; RMS DelaySpread = ' num2str(RMS_DelaySpread*1e9) 'ns ; Rice K = ' num2str(Rice_K)]},'Parent',axes1);
end

return


%% Trial to understand the ISI we face at different BWs in each channel


Filter = 'Gaussian';
%Filter = 'None';

SymbolFs = 2e6;
BT = 0.5;
BW_IF = SymbolFs*8; %4e6 = +-2e6
Duration = 1/SymbolFs*8;

Ts = 1/80e6;
Duration = Nsamples_ChannelAnalysis/80e6;
FTapSize = 80e6/Nsamples_ChannelAnalysis;
FTapSize_MHz = 80/Nsamples_ChannelAnalysis;

%Gaus filter alone:
[Resp_PulseShape, ~, ~, ~] = TimeResp_Gausian(SymbolFs, BT, 1/BW_IF, Duration, FTapSize, BW_IF/2);

%No filter:
if strcmp(Filter ,'None') == 1,
  Resp_PulseShape = [1 zeros(1,length(Resp_PulseShape)-1)];
end


Channel = 3; %0..39

TimeChops = 1:20;%length(Sim_times);
for TimeChop = TimeChops,
    Fcenter = Channel*2 + 2;
    Fstart = Fcenter - BW_IF/1e6/2;
    Fend = Fcenter + BW_IF/1e6/2;
    FirstFTap = Fstart/FTapSize_MHz;
    LastFTap = Fend/FTapSize_MHz;

    ChannelResp = FFChannelResp( TimeChop, FirstFTap+1:LastFTap ); 
    
    figure(1);
    plot(Fstart+FTapSize_MHz:FTapSize_MHz:Fend , 20*log10(abs(ChannelResp./max(ChannelResp))));
    xlabel('MHz + 2400');
    ylabel('dB');
    title([num2str(20*log10(max(abs(ChannelResp)))) 'dBs']);
        
    TimeResp = fftshift(ifft(ChannelResp)); %equivalent to cutting with a rectangular window in freq == like if the Rx filter was a perfect sync
    TimeRespShaped = conv(TimeResp, Resp_PulseShape); %possible pulse shaping
    
    %subsample:
    center = find(TimeRespShaped == max(TimeRespShaped)); %find the maximum
    Step = (BW_IF/SymbolFs); %subsample step
    if mod(center,Step) == 0,
      SampleRange = (-floor(center/Step)+1:(length(TimeRespShaped)-center)/Step)*Step ;
    else
      SampleRange = (-floor(center/Step):(length(TimeRespShaped)-center)/Step)*Step;
    end
    TimeRespShapedSub = TimeRespShaped(center + SampleRange ); %and subsample (Note that subsampling is the same as reducng fs + aliasing)
    
    
    %%thoughts
    [~,center] = max(TimeRespShapedSub);
    ReChannelResp = fft([TimeRespShapedSub(center:end) TimeRespShapedSub(1:center-1)]);
    N= length(ReChannelResp);
    D = sum(ReChannelResp)/N;
    Dd = 20*log10(abs(D));
    AllOther = sqrt(sum(abs(ReChannelResp).^2)/N - abs(D)^2);
    %AllOtherdB = 20*log10(AllOther);
    %fprintf('Freq: DcPower = %2.3f  ; AllOtherPower = %.3f   \n', abs(D)^2, AllOther^2);
    %fprintf('Freq: DcPowdB = %2.3fdB; AllOtherPowdB = %.3fdB \n',  Dd, AllOtherdB);
    %%endofthoughts
    
    %%thoughts
    AllTime = sum(abs(TimeRespShapedSub).^2);
    AllOtherTime = sum(abs(TimeRespShapedSub).^2) - max(abs(TimeRespShapedSub))^2;
    MaxTap= max(abs(TimeRespShapedSub))^2;
    %fprintf('Time: MaxTapPower = %2.3f  ; AllOtherPower = %.3f  ; All   = %.3f \n', MaxTap, AllOtherTime, AllTime);
    %fprintf('Time: MaxTapPowdB = %2.3fdB; AllOtherPowdB = %.3fdB; AlldB = %.3fdB \n\n',10*log10(MaxTap), 10*log10(AllOtherTime), 10*log10(AllTime));
    %%end of thoughts
    
    figure(2);
    hold off;
    Plotcenterd(TimeResp, 1/BW_IF*1e6, 'o', 'DisplayName','TimeResp');
    hold all;
    Plotcenterd(TimeRespShaped, 1/BW_IF*1e6, 'o', 'DisplayName','Shaped');
    Plotcenterd(TimeRespShapedSub, 1/SymbolFs*1e6, 'o', 'DisplayName','Subsampled');
    Plotcenterd(Resp_PulseShape, 1/BW_IF*1e6, 'o', 'DisplayName','ShapFilt');
    legend show
    
    xlabel('us');
    ylabel('dB');
    title([num2str(Sim_times(TimeChop)) 's   ;   ' num2str(20*log10(max(abs(TimeResp)))) 'dB']);
    ylim([-50 0.1]);
    grid on;
    
    %%more thoughts 
    ISI_SNRLoss = 10*log10(MaxTap/AllOtherTime);
    ISI_SNRLossAprox = 20*log10(D/AllOther);
    Fading_LevelLoss = Dd;
    fprintf('Fading_LevelGain = %.3fdB ; ISI = %.3fdB ; ISI_Aprox = %.3fdB \n',Fading_LevelLoss, -ISI_SNRLoss, -ISI_SNRLossAprox);
        
    %%end of more thoughts
    pause
end


return 

