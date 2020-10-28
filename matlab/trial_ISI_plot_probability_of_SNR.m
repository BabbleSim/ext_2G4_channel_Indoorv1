% Copyright 2014 Oticon A/S
% SPDX-License-Identifier: Apache-2.0

opengl software
set(gcf, 'Renderer', 'opengl')

%%
RMS_DelaySpreads = [10, 20, 30 ]*1e-9; %ns
BWs = [2];
Rice_Ks = [0 3]; %For LOS: 0..2 tipical for LOS indoors "small" rooms
            %for NLOS a tipical Rice_K = 0, and shadowing therefore of approx. 3.8dB

            
dBRange = 5:0.5:15.5;
SimTime = 5;%s
DopplerSpeed = 4/3.6; %4km/h
RNGSeed = 1341234;

figure(); hold off;
for Rice_K = Rice_Ks,
  for RMS_DelaySpread = RMS_DelaySpreads,
    for BW = BWs,
      [CDF, bins, pdf] = ISI_SNR_CDF(BW, Rice_K, RMS_DelaySpread, DopplerSpeed, RNGSeed, SimTime, dBRange);


      hplot = plot(bins,CDF*100 ,'DisplayName',[num2str(BW) 'Mb;' num2str(RMS_DelaySpread*1e9) 'ns;K=' num2str(Rice_K) ]);
      %./(BW^2)./(RMS_DelaySpread^2)
      hold all;
      %makedatatip(hplot,find(bins==10,1));
    end
  end
end

legend('Location','Best');
grid on;
xlim([5 max(dBRange)-0.5]);
%ylim([1e-3 0.5]*100);
ylabel('Probability in %');
xlabel('ISI SNR limit in dBs');
title('ISI SNR CDF: probabiltiy of a SNR <= X dBs');
