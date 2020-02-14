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

# This is a generic 128-bit implementation of float to shortest conversion
# using the Ryu algorithm. It can handle any IEEE-compatible floating-point
# type up to 128 bits. In order to use this correctly, you must use the
# appropriate *_to_fd128 function for the underlying type - DO NOT CAST your
# input to another floating-point type, doing so will result in incorrect
# output!
#
# For any floating-point type that is not natively defined by the compiler,
# you can use generic_binary_to_decimal to work directly on the underlying bit
# representation.

const FD128_EXCEPTIONAL_EXPONENT = 0x7FFFFFFF

# A floating decimal representing (-1)^s * m * 10^e.
type floating_decimal_128 = object
  mantissa: uint128
  exponent: int32
  sign: bool

proc float_to_fd128(f: float32): floating_decimal_128
proc double_to_fd128(d: float64): floating_decimal_128

# According to wikipedia (https://en.wikipedia.org/wiki/Long_double), this likely only works on
# x86 with specific compilers (clang?). May need an ifdef.
#proc long_double_to_fd128(d: floatX): floating_decimal_128

# Converts the given binary floating point number to the shortest decimal floating point number
# that still accurately represents it.
proc generic_binary_to_decimal(bits: uint128, mantissaBits, exponentBits: uint32, explicitLeadingBit: bool): floating_decimal_128

# Converts the given decimal floating point number to a string, writing to result, and returning
# the number characters written. Does not terminate the buffer with a 0. In the worst case, this
# function can write up to 53 characters.
#
# Maximal char buffer requirement:
# sign + mantissa digits + decimal dot + 'E' + exponent sign + exponent digits
# = 1 + 39 + 1 + 1 + 1 + 10 = 53
proc generic_to_chars(v: floating_decimal_128, result: var string): int32
