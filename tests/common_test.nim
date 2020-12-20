# Copyright 2019 Ulf Adams
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

suite "common":
  test "decimalLength9":
    check 1u == decimalLength9(0)
    check 1u == decimalLength9(1)
    check 1u == decimalLength9(9)
    check 2u == decimalLength9(10)
    check 2u == decimalLength9(99)
    check 3u == decimalLength9(100)
    check 3u == decimalLength9(999)
    check 9u == decimalLength9(999999999)

  test "ceil_log2pow5":
    check 1 == ceil_log2pow5(0)
    check 3 == ceil_log2pow5(1)
    check 5 == ceil_log2pow5(2)
    check 7 == ceil_log2pow5(3)
    check 10 == ceil_log2pow5(4)
    check 8192 == ceil_log2pow5(3528)

  test "log10Pow2":
    check 0u == log10Pow2(0)
    check 0u == log10Pow2(1)
    check 0u == log10Pow2(2)
    check 0u == log10Pow2(3)
    check 1u == log10Pow2(4)
    check 496u == log10Pow2(1650)

  test "log10Pow5":
    check 0u == log10Pow5(0)
    check 0u == log10Pow5(1)
    check 1u == log10Pow5(2)
    check 2u == log10Pow5(3)
    check 2u == log10Pow5(4)
    check 1831u == log10Pow5(2620)

  test "copy_special_str":
    var buffer: string
    buffer.setLen 3
    check 3 == copy_special_str(buffer, false, false, true)
    check "NaN" == buffer

    buffer.setLen 8
    check 8 == copy_special_str(buffer, false, true, false)
    check "Infinity" == buffer

    buffer.setLen 9
    check 9 == copy_special_str(buffer, true, true, false)
    check "-Infinity" == buffer

    buffer.setLen 3
    check 3 == copy_special_str(buffer, false, false, false)
    check "0E0" == buffer

    buffer.setLen 4
    check 4 == copy_special_str(buffer, true, false, false)
    check "-0E0" == buffer

  test "float_to_bits":
    check 0u == cast[uint32](0.0f)
    check 0x40490fdau == cast[uint32](3.1415926f)

  test "double_to_bits":
    check 0'u64 == cast[uint64](0.0)
    check 0x400921FB54442D18'u64 == cast[uint64](3.1415926535897932384626433)
