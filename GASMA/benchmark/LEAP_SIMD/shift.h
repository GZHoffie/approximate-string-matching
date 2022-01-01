//
// Created by Zhenhao on 2/1/2022.
//

#ifndef LEAP_SIMD_SHIFT_H
#define LEAP_SIMD_SHIFT_H


#ifndef __aligned
#define __aligned __attribute__((aligned(32)))
#endif

#include <x86intrin.h>
#include <stdint.h>


// read and ref need to be 16 aligned
__m128i shift_right_sse(__m128i vec, int shift_num);
__m128i shift_left_sse(__m128i vec, int shift_num);
__m256i shift_right_avx(__m256i vec, int shift_num);
__m256i shift_left_avx(__m256i vec, int shift_num);


#endif //LEAP_SIMD_SHIFT_H
