% Copyright 2014 Oticon A/S
% SPDX-License-Identifier: Apache-2.0

load valuesStore.mat

for trial = 1:length(valuesStore),

Nsamples = size(valuesStore,2);

genvalue = zeros(2^Nsamples,1);

equivalent_std = sqrt(sum( 10.^(valuesStore(trial,:)/10) )) ;

for value = 0:2^Nsamples-1,
  bitmask = 1;
  for bit = 0:Nsamples-1,
    if bitand(bitmask,value),
      genvalue(value+1) = genvalue(value+1) + 10^(valuesStore(trial,bit+1)/20);
    else
      genvalue(value+1) = genvalue(value+1) - 10^(valuesStore(trial,bit+1)/20);
    end
    bitmask = bitmask*2;
  end
end

[pdf, bins] = hist(genvalue,-3*equivalent_std:0.01:3*equivalent_std);
CDF = cumsum(pdf)/length(genvalue);

Equivalent_CDF = alpi_normcdf(bins,equivalent_std);

figure(1);
hold off;
plot(bins,CDF);
hold all
plot(bins,Equivalent_CDF  );

end