function [filter_coefs, FilterCorrection] = Calculate_doppler_filter(Oversampling_recalc)
%
% Copyright 2014 Oticon A/S
% SPDX-License-Identifier: Apache-2.0

%as our sampling period is = Recalc_time = 1/(MaxDopplerSpread*Oversampling_recalc)
%=> fs = (MaxDopplerSpread*Oversampling_recalc) =>
% fs/2 = (MaxDopplerSpread*Oversampling_recalc)/2

%=> the MaxDopplerSpread is placed in 2/Oversampling_recalc (where 1 == fs/2)

if Oversampling_recalc == 16, %let's avoid using licenses. If people doesnt change parameters they dont have to
    filter_coefs = [   0.001949030911796   0.000719809090609   0.000439371138864  -0.000212757947522  -0.001256375437296 ...
                      -0.002641832965118  -0.004239995557018  -0.005845737971487  -0.007199666324979  -0.008021146061509 ...
                      -0.008049136017138  -0.007100346500580  -0.005124630452120  -0.002216635441545   0.001341182269004 ...
                       0.005114281416406   0.008547877799777   0.011042487762756   0.012044583821286   0.011143542178644 ...
                       0.008164960046794   0.003240812441412  -0.003162643620847  -0.010258521384568  -0.017001840803096 ...
                      -0.022207148159751  -0.024687057265213  -0.023419597183090  -0.017701804751496  -0.007278083080951 ...
                       0.007578595327104   0.026053194251744   0.046842312143374   0.068273662539203   0.088482873309502 ...
                       0.105619636111094   0.118060014136404   0.124600883142181   0.124600883142181   0.118060014136404 ...
                       0.105619636111094   0.088482873309502   0.068273662539203   0.046842312143374   0.026053194251744 ...
                       0.007578595327104  -0.007278083080951  -0.017701804751496  -0.023419597183090  -0.024687057265213 ...
                      -0.022207148159751  -0.017001840803096  -0.010258521384568  -0.003162643620847   0.003240812441412 ...
                       0.008164960046794   0.011143542178644   0.012044583821286   0.011042487762756   0.008547877799777 ...
                       0.005114281416406   0.001341182269004  -0.002216635441545  -0.005124630452120  -0.007100346500580 ...
                      -0.008049136017138  -0.008021146061509  -0.007199666324979  -0.005845737971487  -0.004239995557018 ...
                      -0.002641832965118  -0.001256375437296  -0.000212757947522   0.000439371138864   0.000719809090609 ...
                       0.001949030911796 ];
    FilterCorrection = 1 / sqrt(sum(abs(fft(filter_coefs)).^2)/length(filter_coefs));
    filter_coefs = filter_coefs*FilterCorrection;
    FilterCorrection = 1 / sqrt(sum(abs(fft(filter_coefs)).^2)/length(filter_coefs));
else
    Fpass = 2/Oversampling_recalc*0.85;  % Passband Frequency 
    Fstop = 2/Oversampling_recalc*1.25;  % Stopband Frequency
    %== roughly cut at MaxDopplerSpread (doesnt need to be too precise)
    Apass = 1;    % Passband Ripple (dB) == something low enough (we dont care much about the filter complexity)
    Astop = 50;   % Stopband Attenuation (dB) == something big enough
    h = fdesign.lowpass('fp,fst,ap,ast', Fpass, Fstop, Apass, Astop);
    Hd = design(h, 'equiripple');
    filter_coefs =  Hd.Numerator;
    FilterCorrection = 1 / sqrt(sum(abs(fft(filter_coefs)).^2)/length(filter_coefs));
end
