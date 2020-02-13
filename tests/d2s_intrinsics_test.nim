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

import ryu/d2s_intrinsics

suite "d2s_intrinsics":
  test "mod1e9":
    check 0u == mod1e9(0)
    check 1u == mod1e9(1)
    check 2u == mod1e9(2)
    check 10u == mod1e9(10)
    check 100u == mod1e9(100)
    check 1000u == mod1e9(1000)
    check 10000u == mod1e9(10000)
    check 100000u == mod1e9(100000)
    check 1000000u == mod1e9(1000000)
    check 10000000u == mod1e9(10000000)
    check 100000000u == mod1e9(100000000)
    check 0u == mod1e9(1000000000)
    check 0u == mod1e9(2000000000)
    check 1u == mod1e9(1000000001)
    check 1234u == mod1e9(1000001234)
    check 123456789u == mod1e9(12345123456789'u64)
    check 123456789u == mod1e9(123456789123456789'u64)
