% Copyright 2014 Oticon A/S
% SPDX-License-Identifier: Apache-2.0

%% Ideal Rectangular in Freq
BW = 1e6;
Ts = 1/80e6;
Duration = 8e-6;

[Resp, Times] = TimeResp_IdealRect_freq(BW, Ts, Duration);

figure(2);

hold all;
plot(Times*1e6 , 20*log10(abs(Resp./max(Resp))));  ylabel('dB'); ylim([-40 0]);
%plot(Times*1e6 , (abs(Resp./max(Resp))));  ylabel('n.u.'); 
xlabel('us');
grid on;
title('Ideal rect. in freq');


%% Gausian pulse shaping

%http://www.mathworks.se/help/signal/examples/fir-gaussian-pulse-shaping-filter-design.html

BT = 0.5; %B*Ts
SymbolFs = 1e6;

Ts = 1/80e6; %response sampling
Duration = 8e-6;
Fs = 100e3;

[TimeResp, Times, FreqResp, Freqs] = TimeResp_Gausian(SymbolFs, BT, Ts, Duration, Fs, SymbolFs);

plot(Times*1e6,(abs(TimeResp))); xlabel('us');
plot(Freqs/1e6,20*log10(FreqResp)); xlabel('MHz');

%% Rectangle in time
SymbolFs = 1e6;
Ts = 1/80e6;
[Resp_Rect, ~] = TimeResp_Rect(SymbolFs, Ts);


%% 
SymbolFs = 1e6;
BT = 0.5;
BW_IF = 4e6;

Ts = 1/80e6;
Duration = 8e-6;
Fs = 100e3;

[RespIF, ~] = TimeResp_IdealRect_freq(BW_IF, Ts, Duration);
[Resp_Gaus, ~, HGaus, FreqGaus] = TimeResp_Gausian(SymbolFs, BT, Ts, Duration, Fs, Fs);
[Resp_Rect, ~] = TimeResp_Rect(SymbolFs, Ts);

Total = conv(conv(RespIF, Resp_Gaus), Resp_Rect);
figure(1);
hold off;
plot((0:1:length(Total)-1)*Ts*1e6, 20*log10(abs(Total./max(Total))), 'DisplayName','Total'); ylim([-40 0]); xlabel('us'); ylabel('dB');
hold all; 
plot( ( (0:1:length(RespIF)-1) + ( length(Total) -length(RespIF) )/2 ) *Ts*1e6, 20*log10(abs(RespIF./max(RespIF))), 'DisplayName','RespIF'); 
plot( ( (0:1:length(Resp_Gaus)-1) + ( length(Total) -length(Resp_Gaus) )/2 ) *Ts*1e6, 20*log10(abs(Resp_Gaus./max(Resp_Gaus))),'DisplayName','Resp Gaus');
plot( ( (0:1:length(Resp_Rect)-1) + ( length(Total) -length(Resp_Rect) )/2 ) *Ts*1e6, 20*log10(abs(Resp_Rect./max(Resp_Rect))),'DisplayName','Resp Rect');
legend show

figure(2);
hold off;
plot((0:1:length(Total)-1)*Ts*1e6, (abs(Total./max(Total))), 'DisplayName','Total');  xlabel('us'); ylabel('n.u.');
hold all;
plot( ( (0:1:length(RespIF)-1) + ( length(Total) -length(RespIF) )/2 ) *Ts*1e6, (abs(RespIF./max(RespIF))),'DisplayName','RespIF'); 
plot( ( (0:1:length(Resp_Gaus)-1) + ( length(Total) -length(Resp_Gaus) )/2 ) *Ts*1e6, (abs(Resp_Gaus./max(Resp_Gaus))),'DisplayName','Resp Gaus');
plot( ( (0:1:length(Resp_Rect)-1) + ( length(Total) -length(Resp_Rect) )/2 ) *Ts*1e6, (abs(Resp_Rect./max(Resp_Rect))),'DisplayName','Resp Rect');
legend show


%% 

