//
// Created by Zhenhao on 2/1/2022.
//

#ifndef LEAP_SIMD_SHD_H
#define LEAP_SIMD_SHD_H


#ifndef __aligned
#define __aligned __attribute__((aligned(16)))
#endif

#include <stdint.h>
#include <x86intrin.h>

// read and ref need to be 16 aligned
int bit_vec_filter_sse(__m128i read_XMM0, __m128i read_XMM1,
                       __m128i ref_XMM0, __m128i ref_XMM1, int length, int max_error);

int bit_vec_filter_avx(__m256i read_YMM0, __m256i read_YMM1,
                       __m256i ref_YMM0, __m256i ref_YMM1, int length, int max_error);

int bit_vec_filter_avx(__m256i *xor_masks, int length, int max_error);

#endif //LEAP_SIMD_SHD_H
