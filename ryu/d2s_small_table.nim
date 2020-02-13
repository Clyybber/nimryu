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

#NIM:
import common

# Defines HAS_UINT128 and uint128_t if applicable.
import d2s_intrinsics

# These tables are generated by PrintDoubleLookupTable.
const DOUBLE_POW5_INV_BITCOUNT* = 125
const DOUBLE_POW5_BITCOUNT* = 125

const DOUBLE_POW5_INV_SPLIT2: array[13, array[2, uint64]] = [
  [                    1'u64, 2305843009213693952'u64 ],
  [  5955668970331000884'u64, 1784059615882449851'u64 ],
  [  8982663654677661702'u64, 1380349269358112757'u64 ],
  [  7286864317269821294'u64, 2135987035920910082'u64 ],
  [  7005857020398200553'u64, 1652639921975621497'u64 ],
  [ 17965325103354776697'u64, 1278668206209430417'u64 ],
  [  8928596168509315048'u64, 1978643211784836272'u64 ],
  [ 10075671573058298858'u64, 1530901034580419511'u64 ],
  [   597001226353042382'u64, 1184477304306571148'u64 ],
  [  1527430471115325346'u64, 1832889850782397517'u64 ],
  [ 12533209867169019542'u64, 1418129833677084982'u64 ],
  [  5577825024675947042'u64, 2194449627517475473'u64 ],
  [ 11006974540203867551'u64, 1697873161311732311'u64 ]
]
const POW5_INV_OFFSETS: array[19, uint32] = [
  0x54544554'u32, 0x04055545'u32, 0x10041000'u32, 0x00400414'u32, 0x40010000'u32, 0x41155555'u32,
  0x00000454'u32, 0x00010044'u32, 0x40000000'u32, 0x44000041'u32, 0x50454450'u32, 0x55550054'u32,
  0x51655554'u32, 0x40004000'u32, 0x01000001'u32, 0x00010500'u32, 0x51515411'u32, 0x05555554'u32,
  0x00000000'u32
]

const DOUBLE_POW5_SPLIT2: array[13, array[2, uint64]] = [
  [                    0'u64, 1152921504606846976'u64 ],
  [                    0'u64, 1490116119384765625'u64 ],
  [  1032610780636961552'u64, 1925929944387235853'u64 ],
  [  7910200175544436838'u64, 1244603055572228341'u64 ],
  [ 16941905809032713930'u64, 1608611746708759036'u64 ],
  [ 13024893955298202172'u64, 2079081953128979843'u64 ],
  [  6607496772837067824'u64, 1343575221513417750'u64 ],
  [ 17332926989895652603'u64, 1736530273035216783'u64 ],
  [ 13037379183483547984'u64, 2244412773384604712'u64 ],
  [  1605989338741628675'u64, 1450417759929778918'u64 ],
  [  9630225068416591280'u64, 1874621017369538693'u64 ],
  [   665883850346957067'u64, 1211445438634777304'u64 ],
  [ 14931890668723713708'u64, 1565756531257009982'u64 ]
]
const POW5_OFFSETS: array[21, uint32] = [
  0x00000000'u32, 0x00000000'u32, 0x00000000'u32, 0x00000000'u32, 0x40000000'u32, 0x59695995'u32,
  0x55545555'u32, 0x56555515'u32, 0x41150504'u32, 0x40555410'u32, 0x44555145'u32, 0x44504540'u32,
  0x45555550'u32, 0x40004000'u32, 0x96440440'u32, 0x55565565'u32, 0x54454045'u32, 0x40154151'u32,
  0x55559155'u32, 0x51405555'u32, 0x00000105'u32
]

const POW5_TABLE_SIZE = 26
const DOUBLE_POW5_TABLE: array[POW5_TABLE_SIZE, uint64] = [
1'u64, 5'u64, 25'u64, 125'u64, 625'u64, 3125'u64, 15625'u64, 78125'u64, 390625'u64,
1953125'u64, 9765625'u64, 48828125'u64, 244140625'u64, 1220703125'u64, 6103515625'u64,
30517578125'u64, 152587890625'u64, 762939453125'u64, 3814697265625'u64,
19073486328125'u64, 95367431640625'u64, 476837158203125'u64,
2384185791015625'u64, 11920928955078125'u64, 59604644775390625'u64,
298023223876953125'u64 #, 1490116119384765625'u64
]

when defined(HAS_UINT128):
  discard #[ NIM:
  # Computes 5^i in the form required by Ryu, and stores it in the given pointer.
  static inline void double_computePow5(const uint32_t i, uint64_t* const result) {
    const uint32_t base = i / POW5_TABLE_SIZE
    const uint32_t base2 = base * POW5_TABLE_SIZE
    const uint32_t offset = i - base2
    const uint64_t* const mul = DOUBLE_POW5_SPLIT2[base]
    if (offset == 0) {
      result[0] = mul[0]
      result[1] = mul[1]
      return
    }
    const uint64_t m = DOUBLE_POW5_TABLE[offset]
    const uint128_t b0 = ((uint128_t) m) * mul[0]
    const uint128_t b2 = ((uint128_t) m) * mul[1]
    const uint32_t delta = pow5bits(i) - pow5bits(base2)
    const uint128_t shiftedSum = (b0 >> delta) + (b2 << (64 - delta)) + ((POW5_OFFSETS[i / 16] >> ((i % 16) << 1)) & 3)
    result[0] = (uint64_t) shiftedSum
    result[1] = (uint64_t) (shiftedSum >> 64)
  }

  # Computes 5^-i in the form required by Ryu, and stores it in the given pointer.
  static inline void double_computeInvPow5(const uint32_t i, uint64_t* const result) {
    const uint32_t base = (i + POW5_TABLE_SIZE - 1) / POW5_TABLE_SIZE
    const uint32_t base2 = base * POW5_TABLE_SIZE
    const uint32_t offset = base2 - i
    const uint64_t* const mul = DOUBLE_POW5_INV_SPLIT2[base]; # 1/5^base2
    if (offset == 0) {
      result[0] = mul[0]
      result[1] = mul[1]
      return
    }
    const uint64_t m = DOUBLE_POW5_TABLE[offset]; # 5^offset
    const uint128_t b0 = ((uint128_t) m) * (mul[0] - 1)
    const uint128_t b2 = ((uint128_t) m) * mul[1]; # 1/5^base2 * 5^offset = 1/5^(base2-offset) = 1/5^i
    const uint32_t delta = pow5bits(base2) - pow5bits(i)
    const uint128_t shiftedSum =
      ((b0 >> delta) + (b2 << (64 - delta))) + 1 + ((POW5_INV_OFFSETS[i / 16] >> ((i % 16) << 1)) & 3)
    result[0] = (uint64_t) shiftedSum
    result[1] = (uint64_t) (shiftedSum >> 64)
  }
  ]#

else:

  # Computes 5^i in the form required by Ryu, and stores it in the given pointer.
  proc double_computePow5*(i: uint32, result: var array[2, uint64]) {.inline.} =
    let base: uint32 = i div POW5_TABLE_SIZE
    let base2: uint32 = base * POW5_TABLE_SIZE
    let offset: uint32 = i - base2
    let mul: array[2, uint64] = DOUBLE_POW5_SPLIT2[base]
    if offset == 0:
      result[0] = mul[0]
      result[1] = mul[1]
      return
    let m: uint64 = DOUBLE_POW5_TABLE[offset]
    var high1: uint64
    let low1: uint64 = umul128(m, mul[1], high1)
    var high0: uint64
    let low0: uint64 = umul128(m, mul[0], high0)
    let sum: uint64 = high0 + low1
    if sum < high0:
      inc high1; # overflow into high1
    # high1 | sum | low0
    let delta: uint32 = uint32_t pow5bits(int32_t i) - pow5bits(int32_t base2)
    result[0] = shiftright128(low0, sum, delta) + ((POW5_OFFSETS[i div 16] shr ((i mod 16) shl 1)) and 3)
    result[1] = shiftright128(sum, high1, delta)

  # Computes 5^-i in the form required by Ryu, and stores it in the given pointer.
  proc double_computeInvPow5*(i: uint32, result: var array[2, uint64]) {.inline.} =
    let base: uint32 = (i + POW5_TABLE_SIZE - 1) div POW5_TABLE_SIZE
    let base2: uint32 = base * POW5_TABLE_SIZE
    let offset: uint32 = base2 - i
    let mul: array[2, uint64] = DOUBLE_POW5_INV_SPLIT2[base]; # 1/5^base2
    if offset == 0:
      result[0] = mul[0]
      result[1] = mul[1]
      return
    let m: uint64 = DOUBLE_POW5_TABLE[offset]
    var high1: uint64
    let low1: uint64 = umul128(m, mul[1], high1)
    var high0: uint64
    let low0: uint64 = umul128(m, mul[0] - 1, high0)
    let sum: uint64 = high0 + low1
    if sum < high0:
      inc high1; # overflow into high1
    # high1 | sum | low0
    let delta: uint32 = uint32_t pow5bits(int32_t base2) - pow5bits(int32_t i)
    result[0] = shiftright128(low0, sum, delta) + 1 + ((POW5_INV_OFFSETS[i div 16] shr ((i mod 16) shl 1)) and 3)
    result[1] = shiftright128(sum, high1, delta)