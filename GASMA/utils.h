//
// Created by Zhenhao on 4/7/2021.
//

/**
 * This file is adopted from CMU-SAFARI/LEAP
 * https://github.com/CMU-SAFARI/LEAP/blob/master/print.c
 */

#ifndef GASMA_UTILS_H
#define GASMA_UTILS_H

#include <x86intrin.h>
#include <cstdio>

void print_byte_vector(uint8_t *data, int length) {
    int i;
    for (i = 0; i < length; i++) {
        int m;
        for (m = 0; m < 8; m++) {
            if (data[i] & (1ULL << m))
                printf("1");
            else
                printf("0");
        }
    }
}

void print128_bit(__m128i var) {
    uint8_t *val = (uint8_t*) &var;
    print_byte_vector(val, 16);

    printf("\n");
}

void print256_bit(__m256i var) {
    uint8_t *val = (uint8_t*) &var;
    print_byte_vector(val, 32);

    printf("\n");
}


void print128_hex(__m128i var) {
    uint8_t *val = (uint8_t*) &var;
    printf("Numerxcal: %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x\n",
           val[0], val[1], val[2], val[3], val[4], val[5], val[6], val[7],
           val[8], val[9], val[10], val[11], val[12], val[13], val[14],
           val[15]);
}

#endif //GASMA_UTILS_H
