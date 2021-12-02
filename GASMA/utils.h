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
     * Copy constructor of the class `int_128bit` that copies a __m128i object.
     */
    int_128bit(const __m128i & that) {
        val = that;
    }

    /**
     * Copy constructor of the class `int_128bit` that copies a __m128i object.
     */
    int_128bit(const uint8_t * that) {
        val = *((__m128i*) that);
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

    int_128bit _xor(const int_128bit &that) {
        return _mm_xor_si128(this->val, that.val);
    }

    int_128bit _or(const int_128bit &that) {
        return _mm_or_si128(this->val, that.val);
    }

    int_128bit shift_right(int shift_num) {
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

    int_128bit shift_left(int shift_num) {
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

/**
 * Class storing the information of highways.
 */
class highway_info {
public:
    // The starting column of the highway
    int starting_point;

    // Length of the highway
    int length;
};

class hurdle_matrix {
private:
    // maximum difference allowed
    int k;

    // two strings for comparison
    char A[MAX_LENGTH] __aligned__;
    char B[MAX_LENGTH] __aligned__;

    // int_128bit objects storing the bit arrays converted from string A and B
    int_128bit* A_bit0_mask;
    int_128bit* A_bit1_mask;
    int_128bit* B_bit0_mask;
    int_128bit* B_bit1_mask;

    // rows in the hurdle matrix
    int_128bit* lanes;

    // list containing the closest highway of each lane
    class highways {
    private:
        // maximum difference allowed
        int k;
        highway_info* info;

    public:

        /**
         * Constructor of the LSB class.
         * @param error value of parameter k, or the maximum number of difference
         *      allowed.
         */
        explicit highways(int error) {
            k = error;
            info = new highway_info[2 * k + 1];
            info->length = 0;
            info->starting_point = 0;
        }

        highway_info& operator[](int lane) {
            return *info[lane + k];
        }

        ~highways(){
            delete[] info;
        }
    };
    highways* highway_list;

    // current position
    int current_lane;
    int current_column;

    /**
     * Convert A and B into int8 objects and store them in A_bit0_t, A_bit1_t,
     * B_bit0_t and B_bit1_t.
     */
    void _convert_read() {
        // array of int8 objects to store the converted bits
        uint8_t A_bit0_t[MAX_LENGTH / 4] __aligned__;
        uint8_t A_bit1_t[MAX_LENGTH / 4] __aligned__;
        uint8_t B_bit0_t[MAX_LENGTH / 4] __aligned__;
        uint8_t B_bit1_t[MAX_LENGTH / 4] __aligned__;

        // convert string A and B into bits and store in the int8 array
        sse3_convert2bit1(A, A_bit0_t, A_bit1_t);
        sse3_convert2bit1(B, B_bit0_t, B_bit1_t);

        // convert the int8 array into 128-bit integer
        A_bit0_mask = new int_128bit(A_bit0_t);
        A_bit1_mask = new int_128bit(A_bit1_t);
        B_bit0_mask = new int_128bit(B_bit0_t);
        B_bit1_mask = new int_128bit(B_bit1_t);
    }

    void _update_highway_list() {
        for (int lane = -k; lane <= k; lane++) {
            if ((highway_list[lane]).starting_point < current_column) {

            }
        }

    }

    /**
     * Construct 128-bit boolean vectors where the i-th element stores whether
     * read[i] matches ref[i+shift], where shift is the lane id that is between
     * -k and k. Store everything in `lanes`.
     */
    void _construct_hurdles() {
        int_128bit mask_bit0, mask_bit1;
        for (int lane = -k; lane <= k; lane++) {
            if (lane < 0) {
                mask_bit0 = (A_bit0_mask->shift_left(-lane))._xor(*B_bit0_mask);
                mask_bit1 = (A_bit1_mask->shift_left(-lane))._xor(*B_bit1_mask);
            } else {
                mask_bit0 = (B_bit0_mask->shift_left(lane))._xor(*A_bit0_mask);
                mask_bit1 = (B_bit1_mask->shift_left(lane))._xor(*A_bit1_mask);
            }
            auto mask = mask_bit0._or(mask_bit1);
            (*this)[lane] = mask;
        }
    }

public:
    int_128bit& operator[](int lane){
        return lanes[lane + k];
    }

    hurdle_matrix(const char* read, const char* ref, int error, int length) {
        if (length > MAX_LENGTH) {
            length = MAX_LENGTH;
        }
        // assign to class parameters
        strncpy(A, read, length);
        strncpy(B, ref, length);
        k = error;
        _convert_read();
        highway_list = new highways(k);
        lanes = new int_128bit[2 * k + 1];
        _construct_hurdles();

        // define starting position at (0, 0)
        current_lane = 0;
        current_column = 0;
    }

    void print() {
        for (int i = -k; i <= k; i++) {
            (*this)[i].print();
        }
    }

    ~hurdle_matrix() {
        delete highway_list;
        delete[] lanes;
    }
};



#endif //GASMA_UTILS_H
