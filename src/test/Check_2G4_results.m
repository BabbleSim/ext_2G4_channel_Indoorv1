# Copyright 2018 Oticon A/S
# SPDX-License-Identifier: Apache-2.0

opengl software
set(gcf, 'Renderer', 'opengl')

%% Load C model results
clear
load HF.txt
load ChRespF.txt
load AveFade_ISI.txt

params

Modname = Name_from_Modidx(Modulation);
%we load from here Recalc_time & N_h_taps

%% 
Nsamples_ChannelAnalysis = 320; %define in C code
FreqsToPlot = Nsamples_ChannelAnalysis/2; %all or half 
fs = 160e6/ ( Nsamples_ChannelAnalysis/FreqsToPlot ) ;
Ts = 1/160e6;

ChRespFull = ChRespF(:,1) + ChRespF(:,2)*1i;
ChRespFull = reshape(ChRespFull, Nsamples_ChannelAnalysis, []);
%ChResp = fftshift(ChResp ,1);
ChResp = ChRespFull(1:FreqsToPlot,:);
H = HF(:,1) + HF(:,2)*1i;
H = reshape(H, N_h_taps ,[]);
MChRespFull = fft(H ,Nsamples_ChannelAnalysis,1); %we do the FFT ourselves here to compare


Sim_times = (0:size(ChResp,2)-1) * Recalc_time/1e6;
[x,y] = meshgrid( (0:1/FreqsToPlot:1-1/FreqsToPlot)*fs/1e6 , Sim_times);
[x3,y3] = meshgrid( (0:1/Nsamples_ChannelAnalysis:1-1/Nsamples_ChannelAnalysis)*fs/1e6 , Sim_times);

%% Plot the channel impulse response the C model calculated over time
figure(); mesh(max(20*log10(abs(H')),-40)); hold on;
xr = 1:N_h_taps;
yr = [1 size(H,2)];
Values = exp(-Ts.*(xr-1)./DelaySpread);
Values = Values./sqrt(sum(abs(Values).^2));
zs = 20*log10(Values)'*ones(1,2);
[xg,yg] = meshgrid(xr,yr);
mesh(xg,yg,zs','FaceLighting','none',...
  'FaceColor','none',...
  'EdgeColor',[1 0 0]);
figure(); plot(sqrt(sum(abs(H).^2,2)./size(H,2))); hold all; plot(Values); title('mean H from C code vs asymptotic'); %you need to run for a few seconds to converge


%% Plot the channel frequency response the C model calculated 
%plot((0:length(ChResp)-1)/length(ChResp)*160, 20*log10(abs(ChResp)));
ChRespdB = 20*log10(abs(ChResp))';
figure(); mesh(x,y, ChRespdB); hold on ; title('Channel Resp');
figure(); mesh(x3,y3, abs(ChRespFull-MChRespFull)'); hold on ; title('Channel Resp Error');


%% Plot the average fade and ISI SNR the C model calculated
channels = 80;
AveFade = reshape(AveFade_ISI(:,1),channels,[])';
ISI = reshape(AveFade_ISI(:,2),channels,[])';
[x2,y2] = meshgrid( 0:(channels-1) , Sim_times);
figure(); mesh(x2,y2, AveFade); hold on; title(['AveFade, Mod = ' Modname ]);
figure(); mesh(x2,y2, ISI); hold on; title(['ISI, Mod = ' Modname ]);
%figure(); contour(ISI); hold on; title(['ISI, Mod = ' Modname ]);


%% Calculate here the ISI and Average Fading to compare against the C reference and plot the difference


ChRespExtended = [ChRespFull(301:end,:) ; ChRespFull(1:180,:)]'; %we center ourselves with 10 margin around
FTapSize = 0.5;
ChannelCenters = ((0:79) + 10); %we have 10MHz margin on each side
M_ISI = zeros(size(ChRespExtended,1),length(ChannelCenters));
M_AveFade = zeros(size(ChRespExtended,1),length(ChannelCenters));
for ChannelCenter = ChannelCenters,
  Fstarts = ChannelCenter - BW/2;
  Fends = ChannelCenter + BW/2;
  FirstFTap = Fstarts/FTapSize;
  LastFTap = Fends/FTapSize;
  ChannelResp = ChRespExtended( :, FirstFTap+1:LastFTap);
  N= LastFTap - FirstFTap;
  D = abs(sum(ChannelResp,2)/N).^2; %tap 0 power
  AllOther = sum(abs(ChannelResp).^2,2)/N - D; %all others taps powers
  M_ISI(:, ChannelCenter == ChannelCenters) = 10*log10(D./AllOther);
  M_AveFade(:, ChannelCenter == ChannelCenters) = 10*log10(D);
end


figure(); mesh(x2,y2, M_ISI-ISI); hold on; title(['ISI error, Mod = ' Modname ]);
figure(); mesh(x2,y2, M_AveFade-AveFade); hold on; title(['Fade error, Mod = ' Modname ]);

