# Copyright 2018 Oticon A/S
# SPDX-License-Identifier: Apache-2.0

function xo = interpolate(oldtime, x, newtimes)

if length(oldtime) == 1,
  xo  = x * ones(1,length(newtimes));
else
  xo = interp1(oldtime, x, newtimes, 'linear', NaN);
  tmp = interp1(oldtime, x, newtimes, 'nearest', 'extrap');
  xo(isnan(xo)) = tmp(isnan(xo));
end

end
