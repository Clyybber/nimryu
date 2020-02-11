# Copyright 2018 Ulf Adams
#
# The contents of this file may be used under the terms of the Apache License,
# Version 2.0.
#
#    (See accompanying file LICENSE-Apache or copy at
#     http:#www.apache.org/licenses/LICENSE-2.0)
#
# Alternatively, the contents of this file may be used under the terms of
# the Boost Software License, Version 1.0.
#    (See accompanying file LICENSE-Boost or copy at
#     https:#www.boost.org/LICENSE_1_0.txt)
#
# Unless required by applicable law or agreed to in writing, this software
# is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.
#if defined(_M_IX86) || defined(_M_ARM)
#define RYU_32_BIT_PLATFORM
#endif

# Returns the number of decimal digits in v, which must not contain more than 9 digits.
proc decimalLength9(v: uint32): uint32 {.inline.} =
  # Function precondition: v is not a 10-digit number.
  # (f2s: 9 digits are sufficient for round-tripping.)
  # (d2fixed: We print 9-digit blocks.)
  assert v < 1000000000
  if v >= 100000000: return 9
  if v >= 10000000: return 8
  if v >= 1000000: return 7
  if v >= 100000: return 6
  if v >= 10000: return 5
  if v >= 1000: return 4
  if v >= 100: return 3
  if v >= 10: return 2
  return 1

# Returns e == 0 ? 1 : [log_2(5^e)]; requires 0 <= e <= 3528.
proc log2pow5(e: int32): int32 {.inline.} =
  # This approximation works up to the point that the multiplication overflows at e = 3529.
  # If the multiplication were done in 64 bits, it would fail at 5^4004 which is just greater
  # than 2^9297.
  assert e >= 0
  assert e <= 3528
  return (int32_t) ((((uint32_t) e) * 1217359) >> 19);

# Returns e == 0 ? 1 : ceil(log_2(5^e)); requires 0 <= e <= 3528.
proc pow5bits(e: int32): int32 {.inline.} =
  # This approximation works up to the point that the multiplication overflows at e = 3529.
  # If the multiplication were done in 64 bits, it would fail at 5^4004 which is just greater
  # than 2^9297.
  assert e >= 0
  assert e <= 3528
  return (int32_t) (((((uint32_t) e) * 1217359) >> 19) + 1);

# Returns e == 0 ? 1 : ceil(log_2(5^e)); requires 0 <= e <= 3528.
proc ceil_log2pow5(e: int32): int32 {.inline.} =
  return log2pow5(e) + 1

# Returns floor(log_10(2^e)); requires 0 <= e <= 1650.
proc log10Pow2(e: int32): uint32 {.inline.} =
  # The first value this approximation fails for is 2^1651 which is just greater than 10^297.
  assert e >= 0
  assert e <= 1650
  return (((uint32_t) e) * 78913) >> 18;

# Returns floor(log_10(5^e)); requires 0 <= e <= 2620.
proc log10Pow5(e: int32): uint32 {.inline.} =
  # The first value this approximation fails for is 5^2621 which is just greater than 10^1832.
  assert e >= 0
  assert e <= 2620
  return (((uint32_t) e) * 732923) >> 20;

proc copy_special_str(result: var string, sign, exponent, mantissa: bool): int {.inline.} =
  if mantissa:
    memcpy(result, "NaN", 3);
    return 3;
  if sign:
    result[0] = '-';
  if exponent:
    memcpy(result + sign, "Infinity", 8);
    return sign + 8;
  memcpy(result + sign, "0E0", 3);
  return sign + 3;

proc float_to_bits(f: float32): uint32 {.inline.} =
  var bits: uint32 = 0;
  memcpy(addr bits, unsafeAddr f, sizeof(float32));
  return bits;

proc double_to_bits(d: float64): uint64 {.inline.} =
  var bits: uint64 = 0;
  memcpy(addr bits, unsafeAddr d, sizeof(float64));
  return bits;
