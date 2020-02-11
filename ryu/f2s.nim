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

# Runtime compiler options:
# -DRYU_DEBUG Generate verbose debugging output to stdout.

include ryu

include common
include digit_table

when defined(RYU_FLOAT_FULL_TABLE):
  include f2s_full_table
else:
  when defined(RYU_OPTIMIZE_SIZE):
    include d2s_small_table
  else:
    include d2s_full_table
  const FLOAT_POW5_INV_BITCOUNT = (DOUBLE_POW5_INV_BITCOUNT - 64)
  const FLOAT_POW5_BITCOUNT = (DOUBLE_POW5_BITCOUNT - 64)

const FLOAT_MANTISSA_BITS = 23
const FLOAT_EXPONENT_BITS = 8
const FLOAT_BIAS = 127

proc pow5factor_32(value: uint32): uint32 {.inline.} =
  var count: uint32 = 0
  while true:
    assert value != 0
    let q: uint32 = value / 5;
    let r: uint32 = value % 5;
    if r != 0:
      break
    value = q
    inc count
  return count

# Returns true if value is divisible by 5^p.
proc multipleOfPowerOf5_32(value: uint32, p: uint32): bool {.inline.} =
  return pow5factor_32(value) >= p

# Returns true if value is divisible by 2^p.
proc multipleOfPowerOf2_32(value: uint32, p: uint32): bool {.inline.} =
  # __builtin_ctz doesn't appear to be faster here.
  return (value & ((1u << p) - 1)) == 0

# It seems to be slightly faster to avoid uint128_t here, although the
# generated code for uint128_t looks slightly nicer.
proc mulShift32(m: uint32, factor: uint64, shift: int32): uint32 {.inline.} =
  assert shift > 32

  # The casts here help MSVC to avoid calls to the __allmul library
  # function.
  let factorLo: uint32 = (uint32_t)(factor);
  let factorHi: uint32 = (uint32_t)(factor >> 32);
  let bits0: uint64 = (uint64_t)m * factorLo;
  let bits1: uint64 = (uint64_t)m * factorHi;

  when defined(RYU_32_BIT_PLATFORM):
    # On 32-bit platforms we can avoid a 64-bit shift-right since we only
    # need the upper 32 bits of the result and the shift value is > 32.
    let bits0Hi: uint32 = (uint32_t)(bits0 >> 32);
    var bits1Lo: uint32 = (uint32_t)(bits1);
    var bits1Hi: uint32 = (uint32_t)(bits1 >> 32);
    bits1Lo += bits0Hi;
    bits1Hi += (bits1Lo < bits0Hi);
    let s: int32 = shift - 32;
    return (bits1Hi << (32 - s)) | (bits1Lo >> s);
  else: # RYU_32_BIT_PLATFORM
    let sum: uint64 = (bits0 >> 32) + bits1;
    let shiftedSum: uint64 = sum >> (shift - 32);
    assert(shiftedSum <= UINT32_MAX);
    return (uint32_t) shiftedSum;

proc mulPow5InvDivPow2(m: uint32, q: uint32, j: int32): uint32 {.inline.} =
  when defined(RYU_FLOAT_FULL_TABLE):
    return mulShift32(m, FLOAT_POW5_INV_SPLIT[q], j);
  elif defined(RYU_OPTIMIZE_SIZE):
    # The inverse multipliers are defined as [2^x / 5^y] + 1; the upper 64 bits from the double lookup
    # table are the correct bits for [2^x / 5^y], so we have to add 1 here. Note that we rely on the
    # fact that the added 1 that's already stored in the table never overflows into the upper 64 bits.
    var pow5: array[uint64, 2];
    double_computeInvPow5(q, pow5);
    return mulShift32(m, pow5[1] + 1, j);
  else:
    return mulShift32(m, DOUBLE_POW5_INV_SPLIT[q][1] + 1, j);

proc mulPow5divPow2(m: uint32, i: uint32, j: int32): uint32 {.inline.} =
  when defined(RYU_FLOAT_FULL_TABLE):
    return mulShift32(m, FLOAT_POW5_SPLIT[i], j);
  elif defined(RYU_OPTIMIZE_SIZE):
    var pow5: array[uint64, 2];
    double_computePow5(i, pow5);
    return mulShift32(m, pow5[1], j);
  else:
    return mulShift32(m, DOUBLE_POW5_SPLIT[i][1], j);

# A floating decimal representing m * 10^e.
type floating_decimal_32 = object
  mantissa: uint32
  # Decimal exponent's range is -45 to 38
  # inclusive, and can fit in a short if needed.
  exponent: int32

proc f2d(ieeeMantissa: uint32, ieeeExponent: uint32): floating_decimal_32 {.inline.} =
  var e2: int32;
  var m2: uint32;
  if ieeeExponent == 0:
    # We subtract 2 so that the bounds computation has 2 additional bits.
    e2 = 1 - FLOAT_BIAS - FLOAT_MANTISSA_BITS - 2;
    m2 = ieeeMantissa;
  else:
    e2 = (int32_t) ieeeExponent - FLOAT_BIAS - FLOAT_MANTISSA_BITS - 2;
    m2 = (1u << FLOAT_MANTISSA_BITS) | ieeeMantissa;
  let even: bool = (m2 & 1) == 0;
  let acceptBounds: bool = even;

  when defined(RYU_DEBUG):
    printf("-> %u * 2^%d\n", m2, e2 + 2);

  # Step 2: Determine the interval of valid decimal representations.
  let mv: uint32 = 4 * m2;
  let mp: uint32 = 4 * m2 + 2;
  # Implicit bool -> int conversion. True is 1, false is 0.
  let mmShift: uint32 = ieeeMantissa != 0 || ieeeExponent <= 1;
  let mm: uint32 = 4 * m2 - 1 - mmShift;

  # Step 3: Convert to a decimal power base using 64-bit arithmetic.
  var vr, vp, vm: uint32
  var e10: int32
  var vmIsTrailingZeros: bool = false;
  var vrIsTrailingZeros: bool = false;
  var lastRemovedDigit: uint8 = 0;
  if e2 >= 0:
    let q: uint32 = log10Pow2(e2);
    e10 = (int32_t) q;
    let k: int32 = FLOAT_POW5_INV_BITCOUNT + pow5bits((int32_t) q) - 1;
    let i: int32 = -e2 + (int32_t) q + k;
    vr = mulPow5InvDivPow2(mv, q, i);
    vp = mulPow5InvDivPow2(mp, q, i);
    vm = mulPow5InvDivPow2(mm, q, i);
    when defined(RYU_DEBUG):
      printf("%u * 2^%d / 10^%u\n", mv, e2, q);
      printf("V+=%u\nV =%u\nV-=%u\n", vp, vr, vm);
    if (q != 0 && (vp - 1) / 10 <= vm / 10):
      # We need to know one removed digit even if we are not going to loop below. We could use
      # q = X - 1 above, except that would require 33 bits for the result, and we've found that
      # 32-bit arithmetic is faster even on 64-bit machines.
      let l: int32 = FLOAT_POW5_INV_BITCOUNT + pow5bits((int32_t) (q - 1)) - 1;
      lastRemovedDigit = (uint8_t) (mulPow5InvDivPow2(mv, q - 1, -e2 + (int32_t) q - 1 + l) % 10);
    if q <= 9:
      # The largest power of 5 that fits in 24 bits is 5^10, but q <= 9 seems to be safe as well.
      # Only one of mp, mv, and mm can be a multiple of 5, if any.
      if mv % 5 == 0:
        vrIsTrailingZeros = multipleOfPowerOf5_32(mv, q);
      elif acceptBounds:
        vmIsTrailingZeros = multipleOfPowerOf5_32(mm, q);
      else:
        vp -= multipleOfPowerOf5_32(mp, q);
  else:
    let q: uint32 = log10Pow5(-e2);
    e10 = (int32_t) q + e2;
    let i: int32 = -e2 - (int32_t) q;
    let k: int32 = pow5bits(i) - FLOAT_POW5_BITCOUNT;
    var j: int32 = (int32_t) q - k;
    vr = mulPow5divPow2(mv, (uint32_t) i, j);
    vp = mulPow5divPow2(mp, (uint32_t) i, j);
    vm = mulPow5divPow2(mm, (uint32_t) i, j);
    when defined(RYU_DEBUG):
      printf("%u * 5^%d / 10^%u\n", mv, -e2, q);
      printf("%u %d %d %d\n", q, i, k, j);
      printf("V+=%u\nV =%u\nV-=%u\n", vp, vr, vm);
    if (q != 0 && (vp - 1) / 10 <= vm / 10):
      j = (int32_t) q - 1 - (pow5bits(i + 1) - FLOAT_POW5_BITCOUNT);
      lastRemovedDigit = (uint8_t) (mulPow5divPow2(mv, (uint32_t) (i + 1), j) % 10);
    if q <= 1:
      # {vr,vp,vm} is trailing zeros if {mv,mp,mm} has at least q trailing 0 bits.
      # mv = 4 * m2, so it always has at least two trailing 0 bits.
      vrIsTrailingZeros = true;
      if acceptBounds:
        # mm = mv - 1 - mmShift, so it has 1 trailing 0 bit iff mmShift == 1.
        vmIsTrailingZeros = mmShift == 1;
      else:
        # mp = mv + 2, so it always has at least one trailing 0 bit.
        --vp;
    elif q < 31: # TODO(ulfjack): Use a tighter bound here.
      vrIsTrailingZeros = multipleOfPowerOf2_32(mv, q - 1);
      when defined(RYU_DEBUG):
        printf("vr is trailing zeros=%s\n", vrIsTrailingZeros ? "true" : "false");
  when defined(RYU_DEBUG):
    printf("e10=%d\n", e10);
    printf("V+=%u\nV =%u\nV-=%u\n", vp, vr, vm);
    printf("vm is trailing zeros=%s\n", vmIsTrailingZeros ? "true" : "false");
    printf("vr is trailing zeros=%s\n", vrIsTrailingZeros ? "true" : "false");

  # Step 4: Find the shortest decimal representation in the interval of valid representations.
  var removed: int32 = 0;
  var output: uint32
  if vmIsTrailingZeros || vrIsTrailingZeros:
    # General case, which happens rarely (~4.0%).
    while vp / 10 > vm / 10:
      when false: #__clang__ # https:#bugs.llvm.org/show_bug.cgi?id=23106
        # The compiler does not realize that vm % 10 can be computed from vm / 10
        # as vm - (vm / 10) * 10.
        vmIsTrailingZeros &= vm - (vm / 10) * 10 == 0;
      else:
        vmIsTrailingZeros &= vm % 10 == 0;
      vrIsTrailingZeros &= lastRemovedDigit == 0;
      lastRemovedDigit = (uint8_t) (vr % 10);
      vr /= 10;
      vp /= 10;
      vm /= 10;
      inc removed;
    when defined(RYU_DEBUG):
      printf("V+=%u\nV =%u\nV-=%u\n", vp, vr, vm);
      printf("d-10=%s\n", vmIsTrailingZeros ? "true" : "false");
    if vmIsTrailingZeros:
      while vm % 10 == 0:
        vrIsTrailingZeros &= lastRemovedDigit == 0;
        lastRemovedDigit = (uint8_t) (vr % 10);
        vr /= 10;
        vp /= 10;
        vm /= 10;
        inc removed;
    when defined(RYU_DEBUG):
      printf("%u %d\n", vr, lastRemovedDigit);
      printf("vr is trailing zeros=%s\n", vrIsTrailingZeros ? "true" : "false");
    if (vrIsTrailingZeros && lastRemovedDigit == 5 && vr % 2 == 0):
      # Round even if the exact number is .....50..0.
      lastRemovedDigit = 4;
    # We need to take vr + 1 if vr is outside bounds or we need to round up.
    output = vr + ((vr == vm && (!acceptBounds || !vmIsTrailingZeros)) || lastRemovedDigit >= 5);
  else:
    # Specialized for the common case (~96.0%). Percentages below are relative to this.
    # Loop iterations below (approximately):
    # 0: 13.6%, 1: 70.7%, 2: 14.1%, 3: 1.39%, 4: 0.14%, 5+: 0.01%
    while vp / 10 > vm / 10:
      lastRemovedDigit = (uint8_t) (vr % 10);
      vr /= 10;
      vp /= 10;
      vm /= 10;
      inc removed;
    when defined(RYU_DEBUG):
      printf("%u %d\n", vr, lastRemovedDigit);
      printf("vr is trailing zeros=%s\n", vrIsTrailingZeros ? "true" : "false");
    # We need to take vr + 1 if vr is outside bounds or we need to round up.
    output = vr + (vr == vm || lastRemovedDigit >= 5);
  let exp: int32 = e10 + removed;

  when defined(RYU_DEBUG):
    printf("V+=%u\nV =%u\nV-=%u\n", vp, vr, vm);
    printf("O=%u\n", output);
    printf("EXP=%d\n", exp);

  floating_decimal_32 fd;
  fd.exponent = exp;
  fd.mantissa = output;
  return fd;

proc to_chars(v: floating_decimal_32, sign: bool, result: var string): int32 {.inline.} =
  # Step 5: Print the decimal representation.
  var index: int32 = 0;
  if sign:
    result[index] = '-';
    inc index

  var output: uint32 = v.mantissa;
  let olength: uint32 = decimalLength9(output);

  when defined(RYU_DEBUG):
    printf("DIGITS=%u\n", v.mantissa);
    printf("OLEN=%u\n", olength);
    printf("EXP=%u\n", v.exponent + olength);

  # Print the decimal digits.
  # The following code is equivalent to:
  # for (uint32_t i = 0; i < olength - 1; ++i) {
  #   const uint32_t c = output % 10; output /= 10;
  #   result[index + olength - i] = (char) ('0' + c);
  # }
  # result[index] = '0' + output % 10;
  var i: uint32 = 0;
  while output >= 10000:
    when false:#__clang__ # https:#bugs.llvm.org/show_bug.cgi?id=38217
      let c: uint32 = output - 10000 * (output / 10000);
    else:
      let c: uint32 = output % 10000;
    output /= 10000;
    let c0: uint32 = (c % 100) << 1;
    let c1: uint32 = (c / 100) << 1;
    memcpy(result + index + olength - i - 1, DIGIT_TABLE + c0, 2);
    memcpy(result + index + olength - i - 3, DIGIT_TABLE + c1, 2);
    i += 4;
  if output >= 100:
    let c: uint32 = (output % 100) << 1;
    output /= 100;
    memcpy(result + index + olength - i - 1, DIGIT_TABLE + c, 2);
    i += 2;
  if output >= 10:
    let c: uint32 = output << 1;
    # We can't use memcpy here: the decimal dot goes between these two digits.
    result[index + olength - i] = DIGIT_TABLE[c + 1];
    result[index] = DIGIT_TABLE[c];
  else:
    result[index] = (char) ('0' + output);

  # Print decimal point if needed.
  if olength > 1:
    result[index + 1] = '.';
    index += olength + 1;
  else:
    inc index;

  # Print the exponent.
  result[index] = 'E';
  inc index
  var exp: int32 = v.exponent + (int32_t) olength - 1;
  if exp < 0:
    result[index] = '-';
    inc index
    exp = -exp;

  if exp >= 10:
    memcpy(result + index, DIGIT_TABLE + 2 * exp, 2);
    index += 2;
  else:
    result[index] = (char) ('0' + exp);
    inc index

  return index;

proc f2s_buffered_n(f: float32, result: var string): int32 =
  # Step 1: Decode the floating-point number, and unify normalized and subnormal cases.
  const uint32_t bits = float_to_bits(f);

  when defined(RYU_DEBUG):
    # printf("IN=");
    # for (int32_t bit = 31; bit >= 0; --bit) {
    #   printf("%u", (bits >> bit) & 1);
    # }
    # printf("\n");

  # Decode bits into sign, mantissa, and exponent.
  const bool ieeeSign = ((bits >> (FLOAT_MANTISSA_BITS + FLOAT_EXPONENT_BITS)) & 1) != 0;
  const uint32_t ieeeMantissa = bits & ((1u << FLOAT_MANTISSA_BITS) - 1);
  const uint32_t ieeeExponent = (bits >> FLOAT_MANTISSA_BITS) & ((1u << FLOAT_EXPONENT_BITS) - 1);

  # Case distinction; exit early for the easy cases.
  if ieeeExponent == ((1u << FLOAT_EXPONENT_BITS) - 1u) || (ieeeExponent == 0 && ieeeMantissa == 0):
    return copy_special_str(result, ieeeSign, ieeeExponent, ieeeMantissa);

  let v: floating_decimal_32 = f2d(ieeeMantissa, ieeeExponent);
  return to_chars(v, ieeeSign, result);

proc f2s_buffered(f: float32, result: var string) =
  let index: int32 = f2s_buffered_n(f, result);
  # Terminate the string.
  result[index] = '\0';

proc f2s(f: float32): string =
  result.setLen 16
  f2s_buffered(f, result)
