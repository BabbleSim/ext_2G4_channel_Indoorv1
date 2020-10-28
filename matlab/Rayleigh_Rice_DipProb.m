%% Rayleigh calculated probability of fade:
%
% Copyright 2014 Oticon A/S
% SPDX-License-Identifier: Apache-2.0

Nsamples = 1e6;
Rayleigh_Gen = (randn(1,Nsamples) + 1j*randn(1,Nsamples))/sqrt(2);
%hist(abs(Rayleigh_Gen),0:0.02:5);

RaydB = 20*log10(abs(Rayleigh_Gen));
[pdf, bins] = hist(RaydB,-100:0.1:10);
hold off;
%plot(bins,pdf);
CDF = cumsum(pdf)/length(RaydB);
semilogy(-bins,CDF,'DisplayName','calc');
xlim([0 40]);
ylim([1e-4 1]);
grid on;
ylabel('Probability');
xlabel('Fade depth in dBs');
title('Rayleigh (NLOS), probabiltiy of a fade of X dBs or more');


%% theoretical probability of a rayleigh fade
dips = 0:-1:-40; %dBs
Prob = 1-exp(-10.^(dips/10));
hold all;
semilogy(-dips,Prob,'DisplayName','theo');
ylabel('Probability');
xlabel('Fade depth in dBs');
title('Rayleigh (NLOS), probabiltiy of a fade of X dBs or more');
grid on;
legend show;


%% meas probability of fade of a Rice
hold off;
for K = 0:0.5:5
    Nsamples = 1e6;
    
    taps_gain = sqrt( 1 / ( 2* ( 1 + K ) ) ); %== gausians std (not std^2)
    rice_offset = sqrt(K / ( 1 + K)); %in amplitude

    Rice_Gen = (randn(1,Nsamples) + 1j*randn(1,Nsamples))*taps_gain + rice_offset;
    %hist(abs(Rice_Gen),0:0.02:5);

    RicedB = 20*log10(abs(Rice_Gen));
    [pdf, bins] = hist(RicedB,-100:0.1:10);
    %plot(bins,pdf);
    CDF = cumsum(pdf)/length(RicedB);
    semilogy(-bins,CDF,'DisplayName',['K = ' num2str(K)]);
    xlim([0 40]);
    ylim([1e-4 1]);
    grid on;
    hold all;
end
legend show
xlabel('Fade depth in dBs');
ylabel('Probability');
title('Rice, probabiltiy of a fade of X dBs or more');

