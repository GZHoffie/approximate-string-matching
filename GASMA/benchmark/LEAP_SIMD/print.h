//
// Created by Zhenhao on 2/1/2022.
//

#ifndef LEAP_SIMD_PRINT_H
#define LEAP_SIMD_PRINT_H

#include <x86intrin.h>
#include <cstdint>


void printbytevector(uint8_t *data, int length);
void print128_bit(__m128i var);
void print256_bit(__m256i var);
void print128_hex(__m128i var);


#endif //LEAP_SIMD_PRINT_H
