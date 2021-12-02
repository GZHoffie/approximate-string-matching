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


#include <cstdio>
#include <cstring>
#include <x86intrin.h>

#include "bit_convert.h"

// Define constant
#ifndef MAX_LENGTH
#define MAX_LENGTH 128             // The maximum length of string
#endif

void print_byte_vector(uint8_t *data, int length) {
    for (int i = 0; i < length; i++) {
        for (int m = 0; m < 8; m++) {
            if (data[i] & (1ULL << m))
                printf("1");
            else
                printf("0");
        }
    }
}


class int_128bit {
private:
    __m128i val;

public:
    /**
     * Default constructor of the class `int_128bit`. Set everything to 0.
     */
    int_128bit() {
        val = _mm_setr_epi32(0, 0, 0, 0);
    }

    /**
     * Print the value of `val` in binary format.
     */
    void print() {
        auto *val_uint8 = (uint8_t*) &this->val;
        print_byte_vector(val_uint8, 16);
        printf("\n");
    }

    /**
     * Print the value of `val` in hexadecimal format.
     */
    void print_hex() {
        auto *v = (uint8_t*) &this->val;
        for (int i = 0; i < 16; i++) {
            printf("%x", v[i]);
        }
        printf("\n");
    }

    __m128i shift_right(int shift_num) {
        __m128i vec = this->val;
        if (shift_num >= 64) {
            vec = _mm_slli_si128(vec, 8);
            shift_num = shift_num % 64;
        }
        __m128i carryover = _mm_slli_si128(vec, 8);
        carryover = _mm_srli_epi64(carryover, 64 - shift_num);
        vec = _mm_slli_epi64(vec, shift_num);
        return _mm_or_si128(vec, carryover);
    }

    __m128i shift_left(int shift_num) {
        __m128i vec = this->val;
        if (shift_num >= 64) {
            vec = _mm_srli_si128(vec, 8);
            shift_num = shift_num % 64;
        }
        __m128i carryover = _mm_srli_si128(vec, 8);
        carryover = _mm_slli_epi64(carryover, 64 - shift_num);
        vec = _mm_srli_epi64(vec, shift_num);
        return _mm_or_si128(vec, carryover);
    }
};


class hurdle_matrix {
private:
    // maximum difference allowed
    int k;

    // two strings for comparison
    char A[MAX_LENGTH] __aligned__;
    char B[MAX_LENGTH] __aligned__;

    // int8 arrays storing the converted bits
    uint8_t A_bit0_t[MAX_LENGTH / 4] __aligned__;
    uint8_t A_bit1_t[MAX_LENGTH / 4] __aligned__;
    uint8_t B_bit0_t[MAX_LENGTH / 4] __aligned__;
    uint8_t B_bit1_t[MAX_LENGTH / 4] __aligned__;

    // rows in the hurdle matrix
    int_128bit* lanes;

    // list containing the least significant bit of each lane
    class LSB {
    private:
        // maximum difference allowed
        int k;
        unsigned int* position;

    public:
        /**
         * Constructor of the LSB class.
         * @param error value of parameter k, or the maximum number of difference
         *      allowed.
         */
        LSB(int error) {
            k = error;
            position = new unsigned int[2 * k + 1];
        }

        unsigned int& operator[](int lane) {
            return position[lane + k];
        }

        ~LSB(){
            delete[] position;
        }
    };
    LSB* lsb_list;

    /**
     * Convert A and B into int8 objects and store them in A_bit0_t, A_bit1_t,
     * B_bit0_t and B_bit1_t.
     */
    void _convert_read() {
        sse3_convert2bit1(A, A_bit0_t, A_bit1_t);
        sse3_convert2bit1(B, B_bit0_t, B_bit1_t);
    }

    void _update_lsb_list() {

    }

    /**
     * Construct a 128-bit boolean vector where the i-th element stores whether
     * read[i] matches ref[i+shift].
     *
     * @param shift an integer between -k and k.
     * @return boolean vector that stores how read and ref matches.
     */
    int_128bit _construct_hurdles(int shift) {




    }

public:
    int_128bit& operator[](int lane){
        return lanes[lane + k];
    }

    hurdle_matrix(const char* read, const char* ref, unsigned int error, int length) {
        if (length > MAX_LENGTH) {
            length = MAX_LENGTH;
        }
        // assign to class parameters
        strncpy(A, read, length);
        strncpy(B, ref, length);
        k = error;
        _convert_read();
        lsb_list = new LSB(k);

        lanes = new int_128bit[2 * k + 1];
    }

    void print() {
        for (int i = -k; i <= k; i++) {
            (*this)[i].print();
        }
    }

    ~hurdle_matrix() {
        delete lsb_list;
        delete[] lanes;
    }





};



#endif //GASMA_UTILS_H
