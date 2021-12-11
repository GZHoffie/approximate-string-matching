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
#include <cstdlib>
#include <cstring>
#include <x86intrin.h>

#include "bit_convert.h"
#include "mask.h"

// Define constant
#ifndef MAX_LENGTH
#define MAX_LENGTH 128             // The maximum length of string
#endif

#ifndef EXPECTED_ERROR_RATE
#define EXPECTED_ERROR_RATE 1
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
     * Default constructor of the class `int_128bit`. Set value to be 1.
     */
    int_128bit() {
        val = _mm_setr_epi32(1, 0, 0, 0);
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

    /**
     * Perform bit-wise xor with that.
     * @param that an int_128bit object.
     * @return this ^ that
     */
    int_128bit _xor(const int_128bit &that) {
        return _mm_xor_si128(this->val, that.val);
    }

    /**
     * Perform bit-wise or with that.
     * @param that an int_128bit object.
     * @return this | that
     */
    int_128bit _or(const int_128bit &that) {
        return _mm_or_si128(this->val, that.val);
    }

    /**
    * Perform bit-wise and with that.
    * @param that an int_128bit object.
    * @return this & that
    */
    int_128bit _and(const int_128bit &that) {
        return _mm_and_si128(this->val, that.val);
    }

    /**
    * Perform bit-wise not.
    * @return !this
    */
    int_128bit _not() {
        return _mm_andnot_si128(this->val,
                                _mm_setr_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF));
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

    int_128bit shift_right_one() {
        int_128bit one = _mm_setr_epi32(1, 0, 0, 0);
        return this->shift_right(1)._or(one);
    }

    int_128bit shift_left_one() {
        int_128bit reversed_one = _mm_setr_epi32(0, 0, 0, 0x80000000);
        return this->shift_left(1)._or(reversed_one);
    }

    /**
     * Return the index of the lowest set bit.
     */
    int first_one() {
        // TODO: replace this with de Brujin Sequence that is faster than scanning
        uint64_t data [MAX_LENGTH / 64] __aligned__;
        _mm_store_si128((__m128i *) data, this->val);
        int count = 0;
        int trailing_zeros;
        for (uint64_t i : data) {
            trailing_zeros = static_cast<int>(_tzcnt_u64(i));
            count += trailing_zeros;
            if (trailing_zeros < 64) {
                break;
            }
        }
        return count;
    }

    /**
     * Return the index of the lowest unset bit.
     */
    int first_zero() {
        // TODO: replace this with de Brujin Sequence that is faster than scanning
        auto data = this->_not();
        return data.first_one();
    }

    /**
     * Flip the short 1 bit in this->val if both of its neighbors are zeros
     *
     * @return this->val with short 1 bits flipped.
     */
    int_128bit flip_short_hurdles() {
        int_128bit l1 = this->shift_left_one();
        int_128bit r1 = this->shift_right_one();
        int_128bit mask = l1._or(r1);
        return this->_and(mask);
    }

    /**
     * Flip the short 0 bit in this->val if both of its neighbors are ones
     *
     * @return this->val with short 0 bits flipped.
     */
    int_128bit flip_short_matches() {
        int_128bit l1 = this->shift_left_one();
        int_128bit r1 = this->shift_right_one();
        int_128bit mask = l1._and(r1);
        return this->_or(mask);
    }

    /**
     * Count the number of set bits in this->val using hardware POPCNT instruction.
     * Reference: https://stackoverflow.com/questions/17354971/fast-counting-the-number-of-set-bits-in-m128i-register/17355341.
     * @return an integer showing the number of set bits in this->val.
     */
    int pop_count() {
        const __m128i n_hi = _mm_unpackhi_epi64(this->val, this->val);
#ifdef _MSC_VER
        return (int) (__popcnt64(_mm_cvtsi128_si64(n)) + __popcnt64(_mm_cvtsi128_si64(n_hi)));
#else
        return (int) (__popcntq(_mm_cvtsi128_si64(this->val)) + __popcntq(_mm_cvtsi128_si64(n_hi)));
#endif
    }
};


/**
 * Calculate the linear leaping from lane1 to lane2.
 * @param lane1 The lane number that we come from.
 * @param lane2 The lane number that we are going to.
 * @return The leaping penalty.
 */
int linear_leap_lane_penalty(int lane1, int lane2) {
    return abs(lane1 - lane2);
}

/**
 * Calculate the columns skipped while leaping from lane1 to lane2.
 * @param lane1 The lane number that we come from.
 * @param lane2 The lane number that we are going to.
 * @return The number of columns skipped.
 */
int linear_leap_forward_column(int lane1, int lane2) {
    if (lane1 * lane2 >= 0) {
        if (abs(lane1) > abs(lane2)) return abs(lane1) - abs(lane2);
        else return 0;
    }
    return abs(lane1);
}


#endif //GASMA_UTILS_H
