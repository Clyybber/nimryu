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

import std/unittest

import ryu/common
import ryu/d2s_intrinsics
import ryu/d2s_small_table
import ryu/d2s_full_table

suite "d2s_table":
  test "double_computePow5":
    for i in 0'u32..<326'u32:
      var m: array[2, uint64]
      double_computePow5(i, m)
      check m[0] == DOUBLE_POW5_SPLIT[i][0]
      check m[1] == DOUBLE_POW5_SPLIT[i][1]

  test "double_computeInvPow5":
    for i in 0'u32..<292'u32:
      var m: array[2, uint64]
      double_computeInvPow5(i, m)
      check m[0] == DOUBLE_POW5_INV_SPLIT[i][0]
      check m[1] == DOUBLE_POW5_INV_SPLIT[i][1]

