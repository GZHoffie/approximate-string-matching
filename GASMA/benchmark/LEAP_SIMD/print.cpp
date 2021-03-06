//
// Created by Zhenhao on 2/1/2022.
//

#include "print.h"
#include <cstdio>

void printbytevector(uint8_t *data, int length) {
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
    printbytevector(val, 16);

    printf("\n");
}

void print256_bit(__m256i var) {
    uint8_t *val = (uint8_t*) &var;
    printbytevector(val, 32);

    printf("\n");
}


void print128_hex(__m128i var) {
    uint8_t *val = (uint8_t*) &var;
    printf("Numerxcal: %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x\n",
           val[0], val[1], val[2], val[3], val[4], val[5], val[6], val[7],
           val[8], val[9], val[10], val[11], val[12], val[13], val[14],
           val[15]);
}
