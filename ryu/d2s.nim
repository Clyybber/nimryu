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

# Runtime compiler options:
# -DRYU_DEBUG Generate verbose debugging output to stdout.
#
# -DRYU_ONLY_64_BIT_OPS Avoid using uint128_t or 64-bit intrinsics. Slower,
#     depending on your compiler.
#
# -DRYU_OPTIMIZE_SIZE Use smaller lookup tables. Instead of storing every
#     required power of 5, only store every 26th entry, and compute
#     intermediate values with a multiplication. This reduces the lookup table
#     size by about 10x (only one case, and only double) at the cost of some
#     performance. Currently requires MSVC intrinsics.

include ryu

import common
import digit_table
import d2s_intrinsics

#NIM:
import system/ansi_c

# Include either the small or the full lookup tables depending on the mode.
when defined(RYU_OPTIMIZE_SIZE):
  import d2s_small_table
else:
  import d2s_full_table

const DOUBLE_MANTISSA_BITS = 52
const DOUBLE_EXPONENT_BITS = 11
const DOUBLE_BIAS = 1023

proc decimalLength17(v: uint64): uint32 {.inline.} =
  # This is slightly faster than a loop.
  # The average output length is 16.38 digits, so we check high-to-low.
  # Function precondition: v is not an 18, 19, or 20-digit number.
  # (17 digits are sufficient for round-tripping.)
  assert v < 100000000000000000'u64
  if v >= 10000000000000000'u64: return 17
  if v >= 1000000000000000'u64: return 16
  if v >= 100000000000000'u64: return 15
  if v >= 10000000000000'u64: return 14
  if v >= 1000000000000'u64: return 13
  if v >= 100000000000'u64: return 12
  if v >= 10000000000'u64: return 11
  if v >= 1000000000'u64: return 10
  if v >= 100000000'u64: return 9
  if v >= 10000000'u64: return 8
  if v >= 1000000'u64: return 7
  if v >= 100000'u64: return 6
  if v >= 10000'u64: return 5
  if v >= 1000'u64: return 4
  if v >= 100'u64: return 3
  if v >= 10'u64: return 2
  return 1

# A floating decimal representing m * 10^e.
type floating_decimal_64 = object
  mantissa: uint64
  # Decimal exponent's range is -324 to 308
  # inclusive, and can fit in a short if needed.
  exponent: int32

proc d2d(ieeeMantissa: uint64, ieeeExponent: uint32): floating_decimal_64 {.inline.} =
  var e2: int32
  var m2: uint64
  if ieeeExponent == 0:
    # We subtract 2 so that the bounds computation has 2 additional bits.
    e2 = 1 - DOUBLE_BIAS - DOUBLE_MANTISSA_BITS - 2
    m2 = ieeeMantissa
  else:
    e2 = int32_t(ieeeExponent) - DOUBLE_BIAS - DOUBLE_MANTISSA_BITS - 2
    m2 = (1'u64 shl DOUBLE_MANTISSA_BITS) or ieeeMantissa
  let even: bool = (m2 and 1) == 0
  let acceptBounds: bool = even

  when defined(RYU_DEBUG):
    c_printf("-> %llu * 2^%d\n", m2, e2 + 2)

  # Step 2: Determine the interval of valid decimal representations.
  let mv: uint64 = 4 * m2
  # Implicit bool -> int conversion. True is 1, false is 0.
  let mmShift: uint32 = uint32 ord(ieeeMantissa != 0 or ieeeExponent <= 1)
  # We would compute mp and mm like this:
  # uint64_t mp = 4 * m2 + 2
  # uint64_t mm = mv - 1 - mmShift

  # Step 3: Convert to a decimal power base using 128-bit arithmetic.
  var vr, vp, vm: uint64
  var e10: int32
  var vmIsTrailingZeros: bool = false
  var vrIsTrailingZeros: bool = false
  if e2 >= 0:
    # I tried special-casing q == 0, but there was no effect on performance.
    # This expression is slightly faster than max(0, log10Pow2(e2) - 1).
    let q: uint32 = log10Pow2(e2) - uint32 ord(e2 > 3)
    e10 = int32_t(q)
    let k: int32 = DOUBLE_POW5_INV_BITCOUNT + pow5bits(int32_t(q)) - 1
    let i: int32 = -e2 + int32_t(q) + k
    when defined(RYU_OPTIMIZE_SIZE):
      var pow5: array[2, uint64]
      double_computeInvPow5(q, pow5)
      vr = mulShiftAll64(m2, pow5, i, vp, vm, mmShift)
    else:
      vr = mulShiftAll64(m2, DOUBLE_POW5_INV_SPLIT[q], i, vp, vm, mmShift)
    when defined(RYU_DEBUG):
      c_printf("%llu * 2^%d / 10^%u\n", mv, e2, q)
      c_printf("V+=%llu\nV =%llu\nV-=%llu\n", vp, vr, vm)
    if q <= 21:
      # This should use q <= 22, but I think 21 is also safe. Smaller values
      # may still be safe, but it's more difficult to reason about them.
      # Only one of mp, mv, and mm can be a multiple of 5, if any.
      let mvMod5: uint32 = (uint32_t(mv)) - 5 * (uint32_t(div5(mv)))
      if mvMod5 == 0:
        vrIsTrailingZeros = multipleOfPowerOf5(mv, q)
      elif acceptBounds:
        # Same as min(e2 + (~mm & 1), pow5Factor(mm)) >= q
        # <=> e2 + (~mm & 1) >= q && pow5Factor(mm) >= q
        # <=> true && pow5Factor(mm) >= q, since e2 >= q.
        vmIsTrailingZeros = multipleOfPowerOf5(mv - 1 - mmShift, q)
      else:
        # Same as min(e2 + 1, pow5Factor(mp)) >= q.
        vp -= uint64 ord multipleOfPowerOf5(mv + 2, q)
  else:
    # This expression is slightly faster than max(0, log10Pow5(-e2) - 1).
    let q: uint32 = log10Pow5(-e2) - uint32 ord(-e2 > 1)
    e10 = int32_t(q) + e2
    let i: int32 = -e2 - (int32_t) q
    let k: int32 = pow5bits(i) - DOUBLE_POW5_BITCOUNT
    let j: int32 = int32_t(q) - k
    when defined(RYU_OPTIMIZE_SIZE):
      var pow5: array[2, uint64]
      double_computePow5(uint32 i, pow5)
      vr = mulShiftAll64(m2, pow5, j, vp, vm, mmShift)
    else:
      vr = mulShiftAll64(m2, DOUBLE_POW5_SPLIT[i], j, vp, vm, mmShift)
    when defined(RYU_DEBUG):
      c_printf("%llu * 5^%d / 10^%u\n", mv, -e2, q)
      c_printf("%u %d %d %d\n", q, i, k, j)
      c_printf("V+=%llu\nV =%llu\nV-=%llu\n", vp, vr, vm)
    if q <= 1:
      # {vr,vp,vm} is trailing zeros if {mv,mp,mm} has at least q trailing 0 bits.
      # mv = 4 * m2, so it always has at least two trailing 0 bits.
      vrIsTrailingZeros = true
      if acceptBounds:
        # mm = mv - 1 - mmShift, so it has 1 trailing 0 bit iff mmShift == 1.
        vmIsTrailingZeros = mmShift == 1
      else:
        # mp = mv + 2, so it always has at least one trailing 0 bit.
        dec vp
    elif q < 63: # TODO(ulfjack): Use a tighter bound here.
      # We want to know if the full product has at least q trailing zeros.
      # We need to compute min(p2(mv), p5(mv) - e2) >= q
      # <=> p2(mv) >= q && p5(mv) - e2 >= q
      # <=> p2(mv) >= q (because -e2 >= q)
      vrIsTrailingZeros = multipleOfPowerOf2(mv, q)
      when defined(RYU_DEBUG):
        c_printf("vr is trailing zeros=%s\n", if vrIsTrailingZeros: "true" else: "false")
  when defined(RYU_DEBUG):
    c_printf("e10=%d\n", e10)
    c_printf("V+=%llu\nV =%llu\nV-=%llu\n", vp, vr, vm)
    c_printf("vm is trailing zeros=%s\n", if vmIsTrailingZeros: "true" else: "false")
    c_printf("vr is trailing zeros=%s\n", if vrIsTrailingZeros: "true" else: "false")

  # Step 4: Find the shortest decimal representation in the interval of valid representations.
  var removed: int32 = 0
  var lastRemovedDigit: uint8 = 0
  var output: uint64
  # On average, we remove ~2 digits.
  if vmIsTrailingZeros or vrIsTrailingZeros:
    # General case, which happens rarely (~0.7%).
    while true:
      let vpDiv10: uint64 = div10(vp)
      let vmDiv10: uint64 = div10(vm)
      if vpDiv10 <= vmDiv10:
        break
      let vmMod10: uint32 = (uint32_t(vm)) - 10 * (uint32_t(vmDiv10))
      let vrDiv10: uint64 = div10(vr)
      let vrMod10: uint32 = (uint32_t(vr)) - 10 * (uint32_t(vrDiv10))
      vmIsTrailingZeros = vmIsTrailingZeros and vmMod10 == 0
      vrIsTrailingZeros = vrIsTrailingZeros and lastRemovedDigit == 0
      lastRemovedDigit = uint8_t(vrMod10)
      vr = vrDiv10
      vp = vpDiv10
      vm = vmDiv10
      inc removed
    when defined(RYU_DEBUG):
      c_printf("V+=%llu\nV =%llu\nV-=%llu\n", vp, vr, vm)
      c_printf("d-10=%s\n", if vmIsTrailingZeros: "true" else: "false")
    if vmIsTrailingZeros:
      while true:
        let vmDiv10: uint64 = div10(vm)
        let vmMod10: uint32 = (uint32_t(vm)) - 10 * (uint32_t(vmDiv10))
        if vmMod10 != 0:
          break
        let vpDiv10: uint64 = div10(vp)
        let vrDiv10: uint64 = div10(vr)
        let vrMod10: uint32 = (uint32_t(vr)) - 10 * (uint32_t(vrDiv10))
        vrIsTrailingZeros = vrIsTrailingZeros and lastRemovedDigit == 0
        lastRemovedDigit = uint8_t(vrMod10)
        vr = vrDiv10
        vp = vpDiv10
        vm = vmDiv10
        inc removed
    when defined(RYU_DEBUG):
      c_printf("%llu %d\n", vr, lastRemovedDigit)
      c_printf("vr is trailing zeros=%s\n", if vrIsTrailingZeros: "true" else: "false")
    if vrIsTrailingZeros and lastRemovedDigit == 5 and vr mod 2 == 0:
      # Round even if the exact number is .....50..0.
      lastRemovedDigit = 4
    # We need to take vr + 1 if vr is outside bounds or we need to round up.
    output = vr + uint64 ord((vr == vm and (not acceptBounds or not vmIsTrailingZeros)) or lastRemovedDigit >= 5)
  else:
    # Specialized for the common case (~99.3%). Percentages below are relative to this.
    var roundUp: bool = false
    let vpDiv100: uint64 = div100(vp)
    let vmDiv100: uint64 = div100(vm)
    if vpDiv100 > vmDiv100: # Optimization: remove two digits at a time (~86.2%).
      let vrDiv100: uint64 = div100(vr)
      let vrMod100: uint32 = (uint32_t(vr)) - 100 * (uint32_t(vrDiv100))
      roundUp = vrMod100 >= 50
      vr = vrDiv100
      vp = vpDiv100
      vm = vmDiv100
      removed += 2
    # Loop iterations below (approximately), without optimization above:
    # 0: 0.03%, 1: 13.8%, 2: 70.6%, 3: 14.0%, 4: 1.40%, 5: 0.14%, 6+: 0.02%
    # Loop iterations below (approximately), with optimization above:
    # 0: 70.6%, 1: 27.8%, 2: 1.40%, 3: 0.14%, 4+: 0.02%
    while true:
      let vpDiv10: uint64 = div10(vp)
      let vmDiv10: uint64 = div10(vm)
      if vpDiv10 <= vmDiv10:
        break
      let vrDiv10: uint64 = div10(vr)
      let vrMod10: uint32 = (uint32_t(vr)) - 10 * (uint32_t(vrDiv10))
      roundUp = vrMod10 >= 5
      vr = vrDiv10
      vp = vpDiv10
      vm = vmDiv10
      inc removed
    when defined(RYU_DEBUG):
      c_printf("%llu roundUp=%s\n", vr, if roundUp: "true" else: "false")
      c_printf("vr is trailing zeros=%s\n", if vrIsTrailingZeros: "true" else: "false")
    # We need to take vr + 1 if vr is outside bounds or we need to round up.
    output = vr + uint64 ord(vr == vm or roundUp)
  let exp: int32 = e10 + removed

  when defined(RYU_DEBUG):
    c_printf("V+=%llu\nV =%llu\nV-=%llu\n", vp, vr, vm)
    c_printf("O=%llu\n", output)
    c_printf("EXP=%d\n", exp)

  var fd: floating_decimal_64
  fd.exponent = exp
  fd.mantissa = output
  return fd

proc to_chars(v: floating_decimal_64, sign: bool, resul: var string): int {.inline.} =
  # Step 5: Print the decimal representation.
  var index: int32 = 0
  if sign:
    resul[index] = '-'
    inc index

  var output: uint64 = v.mantissa
  let olength: uint32 = decimalLength17(output)

  when defined(RYU_DEBUG):
    c_printf("DIGITS=%llu\n", v.mantissa)
    c_printf("OLEN=%u\n", olength)
    c_printf("EXP=%u\n", uint32_t(v.exponent) + olength)

  # Print the decimal digits.
  # The following code is equivalent to:
  # for (uint32_t i = 0; i < olength - 1; ++i) {
  #   let c: uint32 = output % 10; output /= 10
  #   resul[index + olength - i] = (char) ('0' + c)
  # }
  # resul[index] = '0' + output % 10

  var i: uint32 = 0
  # We prefer 32-bit operations, even on 64-bit platforms.
  # We have at most 17 digits, and uint32_t can store 9 digits.
  # If output doesn't fit into uint32_t, we cut off 8 digits,
  # so the rest will fit into uint32_t.
  if (output shr 32) != 0:
    # Expensive 64-bit division.
    let q: uint64 = div1e8(output)
    var output2: uint32 = (uint32_t(output)) - 100000000 * (uint32_t(q))
    output = q

    let c: uint32 = output2 mod 10000
    output2 = output2 div 10000
    let d: uint32 = output2 mod 10000
    let c0: uint32 = (c mod 100) shl 1
    let c1: uint32 = (c div 100) shl 1
    let d0: uint32 = (d mod 100) shl 1
    let d1: uint32 = (d div 100) shl 1
    resul[(index + int32 olength - i - 1) .. (index + int32 olength - i - 1 + 1)] = cast[string](DIGIT_TABLE[int32(c0)..int32(c0)+1])
    resul[(index + int32 olength - i - 3) .. (index + int32 olength - i - 3 + 1)] = cast[string](DIGIT_TABLE[int32(c1)..int32(c1)+1])
    resul[(index + int32 olength - i - 5) .. (index + int32 olength - i - 5 + 1)] = cast[string](DIGIT_TABLE[int32(d0)..int32(d0)+1])
    resul[(index + int32 olength - i - 7) .. (index + int32 olength - i - 7 + 1)] = cast[string](DIGIT_TABLE[int32(d1)..int32(d1)+1])
    i += 8
  var output2: uint32 = uint32_t(output)
  while output2 >= 10000:
    when false:# __clang__ # https://bugs.llvm.org/show_bug.cgi?id=38217
      let c: uint32 = output2 - 10000 * (output2 / 10000)
    else:
      let c: uint32 = output2 mod 10000
    output2 = output2 div 10000
    let c0: uint32 = (c mod 100) shl 1
    let c1: uint32 = (c div 100) shl 1
    resul[(index + int32 olength - i - 1) .. (index + int32 olength - i - 1 + 1)] = cast[string](DIGIT_TABLE[int32(c0)..int32(c0)+1])
    resul[(index + int32 olength - i - 3) .. (index + int32 olength - i - 3 + 1)] = cast[string](DIGIT_TABLE[int32(c1)..int32(c1)+1])
    i += 4
  if output2 >= 100:
    let c: uint32 = (output2 mod 100) shl 1
    output2 = output2 div 100
    resul[(index + int32 olength - i - 1) .. (index + int32 olength - i - 1 + 1)] = cast[string](DIGIT_TABLE[int32(c)..int32(c)+1])
    i += 2
  if output2 >= 10:
    let c: uint32 = output2 shl 1
    # We can't use memcpy here: the decimal dot goes between these two digits.
    resul[index + int32 olength - i] = DIGIT_TABLE[c + 1]
    resul[index] = DIGIT_TABLE[c]
  else:
    resul[index] = cast[char](uint32_t('0') + output2)

  # Print decimal point if needed.
  if olength > 1:
    resul[index + 1] = '.'
    index += int32 olength + 1
  else:
    inc index

  # Print the exponent.
  resul[index] = 'E'
  inc index
  var exp: int32 = v.exponent + int32_t(olength) - 1
  if exp < 0:
    resul[index] = '-'
    inc index
    exp = -exp

  if exp >= 100:
    let c: int32 = exp mod 10
    resul[index  .. index + 1] = cast[string](DIGIT_TABLE[2 * (exp div 10) .. 2 * (exp div 10) + 1])
    resul[index + 2] = cast[char](int32_t('0') + c)
    index += 3
  elif exp >= 10:
    resul[index  .. index + 1] = cast[string](DIGIT_TABLE[2 * exp .. 2 * exp + 1])
    index += 2
  else:
    resul[index] = cast[char](int32_t('0') + exp)
    inc index

  return index

proc d2d_small_int(ieeeMantissa: uint64, ieeeExponent: uint32, v: var floating_decimal_64): bool {.inline.} =
  let m2: uint64 = (1'u64 shl DOUBLE_MANTISSA_BITS) or ieeeMantissa
  let e2: int32 = int32_t(ieeeExponent) - DOUBLE_BIAS - DOUBLE_MANTISSA_BITS

  if e2 > 0:
    # f = m2 * 2^e2 >= 2^53 is an integer.
    # Ignore this case for now.
    return false

  if e2 < -52:
    # f < 1.
    return false

  # Since 2^52 <= m2 < 2^53 and 0 <= -e2 <= 52: 1 <= f = m2 / 2^-e2 < 2^53.
  # Test if the lower -e2 bits of the significand are 0, i.e. whether the fraction is 0.
  let mask: uint64 = (1'u64 shl -e2) - 1
  let fraction: uint64 = m2 and mask
  if fraction != 0:
    return false

  # f is an integer in the range [1, 2^53).
  # Note: mantissa might contain trailing (decimal) 0's.
  # Note: since 2^53 < 10^16, there is no need to adjust decimalLength17().
  v.mantissa = m2 shr -e2
  v.exponent = 0
  return true

proc d2s_buffered_n(f: float64, resul: var string): int =
  # Step 1: Decode the floating-point number, and unify normalized and subnormal cases.
  let bits: uint64 = double_to_bits(f)

  when defined(RYU_DEBUG):
    discard #NIM:
    # printf("IN=")
    # for (int32_t bit = 63; bit >= 0; --bit) {
    #   printf("%d", (int) ((bits >> bit) & 1))
    # }
    # printf("\n")

  # Decode bits into sign, mantissa, and exponent.
  let ieeeSign: bool = ((bits shr (DOUBLE_MANTISSA_BITS + DOUBLE_EXPONENT_BITS)) and 1) != 0
  let ieeeMantissa: uint64 = bits and ((1'u64 shl DOUBLE_MANTISSA_BITS) - 1)
  let ieeeExponent: uint32 = (uint32_t) ((bits shr DOUBLE_MANTISSA_BITS) and ((1u shl DOUBLE_EXPONENT_BITS) - 1))
  # Case distinction; exit early for the easy cases.
  if ieeeExponent == ((1u shl DOUBLE_EXPONENT_BITS) - 1u) or (ieeeExponent == 0 and ieeeMantissa == 0):
    return copy_special_str(resul, ieeeSign, ieeeExponent != 0, ieeeMantissa != 0)

  var v: floating_decimal_64
  var isSmallInt: bool = d2d_small_int(ieeeMantissa, ieeeExponent, v)
  if isSmallInt:
    # For small integers in the range [1, 2^53), v.mantissa might contain trailing (decimal) zeros.
    # For scientific notation we need to move these zeros into the exponent.
    # (This is not needed for fixed-point notation, so it might be beneficial to trim
    # trailing zeros in to_chars only if needed - once fixed-point notation output is implemented.)
    while true:
      let q: uint64 = div10(v.mantissa)
      let r: uint32 = (uint32_t(v.mantissa)) - 10 * (uint32_t(q))
      if r != 0:
        break
      v.mantissa = q
      inc v.exponent
  else:
    v = d2d(ieeeMantissa, ieeeExponent)

  return to_chars(v, ieeeSign, resul)

proc d2s_buffered(f: float64, result: var string) =
  let index: int = d2s_buffered_n(f, result)

  # Terminate the string.
  #result[index] = '\0'; #NIM: Not for now

#NIM:
from strutils import strip
proc d2s(f: float64): string =
  result.setLen 25
  d2s_buffered(f, result)
  result.strip(chars = {'\0'})
