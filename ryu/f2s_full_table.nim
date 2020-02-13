# Copyright 2018 Ulf Adams
#
# The contents of this file may be used under the terms of the Apache License,
# Version 2.0.
#
#    (See accompanying file LICENSE-Apache or copy at
#     http://www.apache.org/licenses/LICENSE-2.0)
#
# Alternatively, the contents of this file may be used under the terms of
# the Boost Software License, Version 1.0.
#    (See accompanying file LICENSE-Boost or copy at
#     https://www.boost.org/LICENSE_1_0.txt)
#
# Unless required by applicable law or agreed to in writing, this software
# is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.

# This table is generated by PrintFloatLookupTable.
const FLOAT_POW5_INV_BITCOUNT* = 59
const FLOAT_POW5_BITCOUNT* = 61

const FLOAT_POW5_INV_SPLIT*: array[31, uint64] = [
  576460752303423489'u64, 461168601842738791'u64, 368934881474191033'u64, 295147905179352826'u64,
  472236648286964522'u64, 377789318629571618'u64, 302231454903657294'u64, 483570327845851670'u64,
  386856262276681336'u64, 309485009821345069'u64, 495176015714152110'u64, 396140812571321688'u64,
  316912650057057351'u64, 507060240091291761'u64, 405648192073033409'u64, 324518553658426727'u64,
  519229685853482763'u64, 415383748682786211'u64, 332306998946228969'u64, 531691198313966350'u64,
  425352958651173080'u64, 340282366920938464'u64, 544451787073501542'u64, 435561429658801234'u64,
  348449143727040987'u64, 557518629963265579'u64, 446014903970612463'u64, 356811923176489971'u64,
  570899077082383953'u64, 456719261665907162'u64, 365375409332725730'u64
]
const FLOAT_POW5_SPLIT*: array[47, uint64] = [
  1152921504606846976'u64, 1441151880758558720'u64, 1801439850948198400'u64, 2251799813685248000'u64,
  1407374883553280000'u64, 1759218604441600000'u64, 2199023255552000000'u64, 1374389534720000000'u64,
  1717986918400000000'u64, 2147483648000000000'u64, 1342177280000000000'u64, 1677721600000000000'u64,
  2097152000000000000'u64, 1310720000000000000'u64, 1638400000000000000'u64, 2048000000000000000'u64,
  1280000000000000000'u64, 1600000000000000000'u64, 2000000000000000000'u64, 1250000000000000000'u64,
  1562500000000000000'u64, 1953125000000000000'u64, 1220703125000000000'u64, 1525878906250000000'u64,
  1907348632812500000'u64, 1192092895507812500'u64, 1490116119384765625'u64, 1862645149230957031'u64,
  1164153218269348144'u64, 1455191522836685180'u64, 1818989403545856475'u64, 2273736754432320594'u64,
  1421085471520200371'u64, 1776356839400250464'u64, 2220446049250313080'u64, 1387778780781445675'u64,
  1734723475976807094'u64, 2168404344971008868'u64, 1355252715606880542'u64, 1694065894508600678'u64,
  2117582368135750847'u64, 1323488980084844279'u64, 1654361225106055349'u64, 2067951531382569187'u64,
  1292469707114105741'u64, 1615587133892632177'u64, 2019483917365790221'u64
]

