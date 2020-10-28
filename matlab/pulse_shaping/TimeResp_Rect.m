function [TimeResp, Times] = TimeResp_Rect(SymbolFs, Ts)
%
% Copyright 2014 Oticon A/S
% SPDX-License-Identifier: Apache-2.0

TimeResp = ones(1,(1/SymbolFs)/Ts);
Times = (0:1:(1/SymbolFs)/Ts-1)*Ts;
