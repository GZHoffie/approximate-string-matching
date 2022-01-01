//
// Created by Zhenhao on 2/1/2022.
//

#ifndef LEAP_SIMD_BIT_CONVERT_H
#define LEAP_SIMD_BIT_CONVERT_H


#include <stdint.h>

#ifndef __aligned
#define __aligned __attribute__((aligned(32)))
#endif

void c_convert2bit(char *str, int length, uint8_t *bits);

void sse_convert2bit(char *str, uint8_t *bits0, uint8_t *bits1);

void avx_convert2bit(char *str, uint8_t *bits0, uint8_t *bits1);


#endif //LEAP_SIMD_BIT_CONVERT_H
