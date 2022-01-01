//
// Created by Zhenhao on 2/1/2022.
//

#ifndef LEAP_SIMD_POPCOUNT_H
#define LEAP_SIMD_POPCOUNT_H

#ifndef __aligned
#define __aligned __attribute__((aligned(32)))
#endif

#include <x86intrin.h>
#include <cstdint>


uint32_t popcount_m128i_sse(__m128i reg);
uint32_t popcount_m256i_avx(__m256i reg);

// SHD popcount with SRS mask.
uint32_t popcount_SHD_sse(__m128i reg);
uint32_t popcount_SHD_avx(__m256i reg);

uint32_t builtin_popcount(uint8_t* buffer, int chunks16);

uint32_t popcount(uint8_t *buffer, int chunks16);


#endif //LEAP_SIMD_POPCOUNT_H
