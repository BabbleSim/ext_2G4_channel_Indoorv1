%
% Copyright 2014 Oticon A/S
% SPDX-License-Identifier: Apache-2.0

S = 32e-9;
t= 0:1e-9:S*10;
h = exp(-t/S);
plot(t*1e9+30,20*log10(h))
ylim([-60 20])
xlim([0 200]);
ylabel('db');
xlabel('ns');
title({'Average impulse response with S =32ns','(delayed 30ns for comparison)'});
grid on


%%
tau=[0 50 110 170 290 310]*1e-9;
        % vector of average path power gains in dB
pdb=[0 -3 -10 -18 -26 -32];
pdb = pdb - 10*log10(sum(10.^(pdb./10)));
hold all; plot(tau,pdb);