function [values] = alpi_normcdf(x, std)
%
% To not pay for an extra license just for a 1 liner function
%
% Copyright 2014 Oticon A/S
% SPDX-License-Identifier: Apache-2.0

values = 1/2*(1+erf(x/sqrt(2)/std));