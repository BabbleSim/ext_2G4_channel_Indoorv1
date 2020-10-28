% Copyright 2014 Oticon A/S
% SPDX-License-Identifier: Apache-2.0

model = 'ITU';
Oversampling = 4; %relative to 1 channel
Ts = 1e-6/Oversampling;
RmsDelay = 150e-9; %indoor residential with delay spread in the range of 20 to 150ns
                  %big open areas up to 500ns in extreme cases
            %rms delay corresponds to a e^(-1) decay in the power delay
            %profile == -8.68dB. A 30dB decay == takes 3.45 times longer
               
DopplerSpread = 10; %below 10Hz for walking granpas: 2.4e9 * 1.1m/s (4km/h) / c

switch model
    case 'ITU',
        %  vector of path delays in sec
        tau=[0 50 110 170 290 310]*1e-9;
        % vector of average path power gains in dB
        pdb=[0 -3 -10 -18 -26 -32] - 10*log10(sum(10.^(pdb./10)));
        k = 0;
    case 'alpi',
        %  vector of path delays in sec
        tau=[0 50 110 170 290 310]*1e-9 / (300e-9/(-30/(20*log10(exp(-1)))*RmsDelay)) ;
        % vector of average path power gains in dB
        pdb=[0 -3 -10 -18 -26 -32];
        pdb = pdb - 10*log10(sum(10.^(pdb./10)));
        %plot(tau*1e9,pdb)
        k = 0; %k = 0 == No LOS
               %if loss, seems around k=1.4 is reasonable if it was only 1 tap
    otherwise
        error('Model not implemented.');
end

% Object for Rayleigh channel with classic Doppler spectrum
h = ricianchan(Ts,DopplerSpread,k,tau,pdb);

h.DopplerSpectrum=doppler.flat;
h.ResetBeforeFiltering=0;
h.StoreHistory = true;

tx = randint(500, 1, 2);

dpskSig = fskmod(tx, 2, 500e3 , Oversampling, 1/Ts);
y = filter(h, dpskSig);
plot(h);    


%%

model = 'alpi';
Oversampling = 4*20; %relative to 1 channel
Ts = 1e-6/Oversampling;
RmsDelay = 70e-9; %indoor residential with delay spread in the range of 20 to 150ns
                  %big open areas up to 500ns in extreme cases
            %rms delay corresponds to a e^(-1) decay in the power delay
            %profile == -8.68dB. A 30dB decay == takes 3.45 times longer
            %86 ns == ITU default
               
DopplerSpread = 10e3; %below 10Hz for walking granpas: 2.4e9 * 1.1m/s (4km/h) / c

switch model
    case 'ITU',
        %  vector of path delays in sec
        tau=[0 50 110 170 290 310]*1e-9;
        % vector of average path power gains in dB
        pdb=[0 -3 -10 -18 -26 -32] - 10*log10(sum(10.^(pdb./10)));
        k = 0;
    case 'alpi',
        %  vector of path delays in sec
        tau=[0 50 110 170 290 310]*1e-9 / (300e-9/(-30/(20*log10(exp(-1)))*RmsDelay)) ;
        % vector of average path power gains in dB
        pdb=[0 -3 -10 -18 -26 -32];
        pdb = pdb - 10*log10(sum(10.^(pdb./10)));
        %plot(tau*1e9,pdb)
        k = 0; %k = 0 == No LOS
               %if loss, seems around k=1.4 is reasonable if it was only 1 tap
    otherwise
        error('Model not implemented.');
end

% Object for Rayleigh channel with classic Doppler spectrum
h = ricianchan(Ts,DopplerSpread,k,tau,pdb);

h.DopplerSpectrum=doppler.flat;
h.ResetBeforeFiltering=0;
h.StoreHistory = true;

duration = 1e-4;
Nsamples = duration/Ts;
inputsignal = ones(1,Nsamples); %anything to pass thru
y = filter(h, inputsignal);
plot(h);