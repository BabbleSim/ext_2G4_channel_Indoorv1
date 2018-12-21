# Copyright 2018 Oticon A/S
# SPDX-License-Identifier: Apache-2.0

function  [name] = Name_from_Modidx(idx)

switch idx,
  case 1
    name = 'BLE';
  case 2
    name = 'Prop2Mv1';
  case 3
    name = 'Prop3Mv1';
  case 4
    name = 'Prop4Mv1';
  case 256
    name = 'WLANInter';
  case 257
    name = 'BTInter';
  case 258
    name = 'BLEInter';
  case 259
    name = 'CWInter';
end