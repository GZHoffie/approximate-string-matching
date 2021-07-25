//
// Created by Zhenhao on 4/7/2021.
//

/**
 * This file is adopted from CMU-SAFARI/LEAP
 * https://github.com/CMU-SAFARI/LEAP/blob/master/print.c
 * https://github.com/CMU-SAFARI/LEAP/blob/master/shift.c
 */

#ifndef GASMA_UTILS_H
#define GASMA_UTILS_H

#include <x86intrin.h>
#include <cstdio>

// Define constant
#define _MAX_LENGTH_ 128             // The maximum length of string






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


class int_128bit {
private:
    __m128i val;

public:
    /**
     * Constructor of the class `int_128bit`.
     *
     * @param read the DNA string consisting of `A`, `C`, `D` or `T`,
     *      and convert
     */
    int_128bit(std::string read) {

    }

    /**
     * Print the value of `val` in binary format.
     */
    void print() {
        uint8_t *val_uint8 = (uint8_t*) &this->val;
        print_byte_vector(val_uint8, 16);
        printf("\n");
    }

    /**
     * Print the value of `val` in hexadecimal format.
     */
    void print_hex() {
        uint8_t *val = (uint8_t*) &this->val;
        printf("%x%x%x%x %x%x%x%x %x%x%x%x %x%x%x%x\n",
               val[0], val[1], val[2], val[3], val[4], val[5], val[6], val[7],
               val[8], val[9], val[10], val[11], val[12], val[13], val[14],
               val[15]);
    }

    __m128i shift_right(int shift_num) {
        if (shift_num >= 64) {
            vec = _mm_slli_si128(this->val, 8);
            shift_num = shift_num % 64;
        }
        __m128i carryover = _mm_slli_si128(this->val, 8);
        carryover = _mm_srli_epi64(carryover, 64 - shift_num);
        vec = _mm_slli_epi64(this->val, shift_num);
        return _mm_or_si128(this->val, carryover);
    }

    __m128i shift_left(int shift_num) {
        if (shift_num >= 64) {
            vec = _mm_srli_si128(this->val, 8);
            shift_num = shift_num % 64;
        }
        __m128i carryover = _mm_srli_si128(this->val, 8);
        carryover = _mm_slli_epi64(carryover, 64 - shift_num);
        vec = _mm_srli_epi64(this->val, shift_num);
        return _mm_or_si128(this->val, carryover);
    }
};

#endif //GASMA_UTILS_H
