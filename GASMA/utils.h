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
     * Flip the short 1 bit in this->val if
     * 1. both of its neighbors are zeros
     * 2. both of its 2nd neighbors are zeros or
     *    either of its 2nd neighbors and either of its 3rd neighbors are zeros.
     *
     * @return this->val with short 1 bits flipped.
     */
    int_128bit flip_short_hurdles() {
        int_128bit l1 = this->shift_left_one();
        int_128bit l2 = l1.shift_left_one();
        int_128bit l3 = l2.shift_left_one();
        int_128bit r1 = this->shift_right_one();
        int_128bit r2 = r1.shift_right_one();
        int_128bit r3 = l2.shift_right_one();

        int_128bit both_neighbor_zero = l1._or(r1);
        int_128bit both_2nd_neighbor_zero = l2._or(r2);
        int_128bit either_2nd_neighbor_zero = l2._and(r2);
        int_128bit either_3rd_neighbor_zero = l3._and(r3);

        int_128bit mask_1 = both_neighbor_zero._or(both_2nd_neighbor_zero);
        int_128bit mask_2 = both_neighbor_zero._or(either_2nd_neighbor_zero)._or(either_3rd_neighbor_zero);


        return this->_and(mask_1._and(mask_2));
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
 * Class storing the information of highways.
 */
class highway_info {
public:
    // The starting column of the highway
    int starting_point;

    // Length of the highway
    int length;

    // Cost to reach this highway
    int cost;

    // Destination column
    int destination;

    // Position of first hurdle
    int first_hurdle;
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


class hurdle_matrix {
private:
    // maximum difference allowed
    int k;

    // two strings for comparison
    char A[MAX_LENGTH] __aligned__;
    char B[MAX_LENGTH] __aligned__;

#ifdef DEBUG
    // string storing the original two strings
    char A_orig[MAX_LENGTH] __aligned__;
    char B_orig[MAX_LENGTH] __aligned__;

    // strings storing the matched strings
    char A_match[MAX_LENGTH * 2] __aligned__;
    char B_match[MAX_LENGTH * 2] __aligned__;
    int A_index, B_index, A_match_index, B_match_index;
#endif

    // int_128bit objects storing the bit arrays converted from string A and B
    int_128bit* A_bit0_mask;
    int_128bit* A_bit1_mask;
    int_128bit* B_bit0_mask;
    int_128bit* B_bit1_mask;

    // rows in the hurdle matrix
    int_128bit* lanes;
    int_128bit* lanes_without_short_hurdles;

    // information about destination
    int destination_lane;

    // length of read and ref
    int m, n;

    // CIGAR string used in mapper
    char CIGAR[MAX_LENGTH * 2];
    int CIGAR_index;

#ifdef DEBUG
    /**
     * Update the matching strings given the lanes we are leaping to and the distance we
     * are moving.
     * @param best_lane The lane we are leaping to.
     * @param curr_lane The lane we are currently at.
     * @param distance The distance we travel on best_lane to reach the highway.
     */
    void _update_match(int best_lane, int curr_lane, int distance) {
        // update matched strings
        if (best_lane < curr_lane) {
            for (int i = best_lane; i < curr_lane; i++) {
                A_match[A_match_index] = A_orig[A_index];
                B_match[B_match_index] = '-';
                A_match_index++, B_match_index++, A_index++;
            }
        } else {
            for (int i = curr_lane; i < best_lane; i++) {
                A_match[A_match_index] = '-';
                B_match[B_match_index] = B_orig[B_index];
                A_match_index++, B_match_index++, B_index++;
            }
        }
        for (int i = 0; i < distance; i++) {
            A_match[A_match_index] = A_orig[A_index];
            B_match[B_match_index] = B_orig[B_index];
            A_match_index++, B_match_index++, A_index++, B_index++;
        }
        printf("%s\n%s\n", A_match, B_match);
    }
#endif

    /**
     * Update the CIGAR strings given the lanes we are leaping to and the distance we
     * are moving.
     * @param best_lane The lane we are leaping to.
     * @param curr_lane The lane we are currently at.
     * @param distance The distance we travel on best_lane to reach the highway.
     */
    void _update_CIGAR(int best_lane, int curr_lane, int distance) {
        // update matched strings
        if (best_lane < curr_lane) {
            for (int i = best_lane; i < curr_lane; i++) {
                CIGAR[CIGAR_index] = 'D';
                CIGAR_index++;
            }
        } else {
            for (int i = curr_lane; i < best_lane; i++) {
                CIGAR[CIGAR_index] = 'I';
                CIGAR_index++;
            }
        }
        for (int i = 0; i < distance; i++) {
            CIGAR[CIGAR_index] = 'M';
            CIGAR_index++;
        }
    }

    // list containing the closest highway of each lane
    class highways {
    private:
        // maximum difference allowed
        int k;
        highway_info* info;

        static int _calculate_destination(int _m, int _n, int lane) {
            if (_m >= _n) {
                if (lane < 0) return _n + lane;
                else if (lane < _m - _n) return _n;
                else return _m - lane;
            } else {
                if (lane > 0) return _m - lane;
                else if (lane < _m - _n) return _m;
                else return _n - lane;
            }
        }

    public:
        int best_highway_lane;

        /**
         * Constructor of the LSB class.
         * @param error value of parameter k, or the maximum number of difference
         *      allowed.
         * @param m, n the length of read and ref string
         */
        explicit highways(int error, int m, int n) {
            k = error;
            info = new highway_info[2 * k + 1];
            best_highway_lane = 0;
            for (int lane = -k; lane <= k; lane++) {
                info[lane + k].starting_point = -1;
                info[lane + k].cost = MAX_LENGTH;
                info[lane + k].destination = _calculate_destination(m, n, lane) + 1;
                info[lane + k].first_hurdle = -1;
            }
        }

        highway_info& operator[](int lane) {
            return info[lane + k];
        }

        void print() {
            for (int lane = -k; lane <= k; lane++) {
                printf("lane: %d, starting point: %d, length: %d, cost: %d\n", lane,
                       (*this)[lane].starting_point, (*this)[lane].length, (*this)[lane].cost);
            }
        }

        ~highways(){
            delete[] info;
        }
    };
    highways* highway_list;

    // current position
    int current_lane;
    int current_column;

    // total cost
    int cost;

    /**
     * Convert A and B into int8 objects and store them in A_bit0_t, A_bit1_t,
     * B_bit0_t and B_bit1_t.
     */
    void _convert_read() {
        // array of int8 objects to store the converted bits
        uint8_t A_bit0_t[MAX_LENGTH / 8] __aligned__;
        uint8_t A_bit1_t[MAX_LENGTH / 8] __aligned__;
        uint8_t B_bit0_t[MAX_LENGTH / 8] __aligned__;
        uint8_t B_bit1_t[MAX_LENGTH / 8] __aligned__;

        // convert string A and B into bits and store in the int8 array
        sse3_convert2bit1(A, A_bit0_t, A_bit1_t);
        sse3_convert2bit1(B, B_bit0_t, B_bit1_t);

        // convert the int8 array into 128-bit integer
        A_bit0_mask = new int_128bit(A_bit0_t);
        A_bit1_mask = new int_128bit(A_bit1_t);
        B_bit0_mask = new int_128bit(B_bit0_t);
        B_bit1_mask = new int_128bit(B_bit1_t);
    }

    /**
     * Count the effective length
     * @param lane
     * @param column
     * @return
     */
    int _count_effective_length(int lane, int column) {
        // Extract the highway
        return 0;

    }

    /**
     * Update highway_list for each lane and show the closest highway to the
     * current position. Calculate the cost and record the longest highway.
     *
     * @return a boolean value indicating whether there is still a valid highway.
     */
    bool _update_highway_list() {
        int smallest_highway_heuristic = INT32_MAX;
        int best_highway_lane = 0;
        int first_zero;
        bool reaching_destination = false; // check if we are reaching destination
        for (int lane = -k; lane <= k; lane++) {
            int start_col = current_column + linear_leap_forward_column(current_lane, lane);
            if ((*highway_list)[lane].starting_point < start_col) {
                // get closest highway in the lane
                int_128bit l = (lanes[lane + k]).shift_left(start_col);
                int_128bit l_without_hurdles = (lanes_without_short_hurdles[lane + k]).shift_left(start_col);

                // update highway in lane
                first_zero = l.first_zero();
                int next_hurdle = (l.shift_left(first_zero)).first_one();
                int length = (l_without_hurdles.shift_left(first_zero)).first_one();
                (*highway_list)[lane].starting_point = start_col + first_zero;
                (*highway_list)[lane].length = length;
                (*highway_list)[lane].first_hurdle = next_hurdle;

                // Fix length if reaches destination
                if (start_col + first_zero + length > (*highway_list)[lane].destination) {
                    (*highway_list)[lane].length = std::max(0, (*highway_list)[lane].destination -\
                                                   (start_col + first_zero) + 1);
                    (*highway_list)[lane].first_hurdle = std::min((*highway_list)[lane].first_hurdle,
                                                                  (*highway_list)[lane].length);
                    reaching_destination = true;
                }
            }
            // calculate cost to reach the highway
            // FIXME: use function pointer for more complicated penalties.
            int leap_cost = linear_leap_lane_penalty(current_lane, lane);
            int hurdle_cost = (*highway_list)[lane].starting_point - start_col;
            (*highway_list)[lane].cost = leap_cost + hurdle_cost;

            // calculate expected cost to reach destination
            //leap_cost = linear_leap_lane_penalty(lane, destination_lane);
            //hurdle_cost = (*highway_list)[lane].destination - (*highway_list)[lane].starting_point -
            //        (*highway_list)[lane].length;
            //(*highway_list)[lane].heuristic = leap_cost + (int)(EXPECTED_ERROR_RATE * hurdle_cost);

        }
        for (int lane = -k; lane <= k; lane++) {
            // get the best-looking highway
            int heuristic = (*highway_list)[lane].cost - (*highway_list)[lane].length;
            if (reaching_destination) {
                heuristic += linear_leap_lane_penalty(lane, destination_lane);
            }
            if (heuristic < smallest_highway_heuristic) {
                smallest_highway_heuristic = heuristic;
                best_highway_lane = lane;
            }
        }
        highway_list->best_highway_lane = best_highway_lane;
#ifdef DEBUG
        highway_list->print();
#endif
        if ((*highway_list)[best_highway_lane].length <= 0) {
            return false;
        }
        return true;
    }

    /**
     * Choose the best highway according to the updated highway list.
     * @return The lane number where the best highway is.
     */
    int _choose_best_highway() {
        int best_lane = highway_list->best_highway_lane;
        int starting_point = (*highway_list)[best_lane].starting_point;
        int smallest_cost = (*highway_list)[best_lane].cost;
        int smallest_cost_length = 0;
        int smallest_cost_lane = best_lane;
        for (int lane = -k; lane <= k; lane++) {
            if (lane != best_lane) {
                int ending_point = (*highway_list)[lane].starting_point + (*highway_list)[lane].length;
                int new_cost = (*highway_list)[lane].cost + linear_leap_lane_penalty(lane, best_lane) \
                               + std::max(0, starting_point - linear_leap_forward_column(lane, best_lane) - ending_point);
                if (new_cost < smallest_cost || (new_cost == smallest_cost &&
                    (*highway_list)[lane].length > smallest_cost_length)) {
                    smallest_cost = new_cost;
                    smallest_cost_lane = lane;
                    smallest_cost_length = (*highway_list)[lane].length;
                }
            }
        }
        return smallest_cost_lane;
    }

    /**
     * Perform one step in the greedy algorithm.
     * @return a boolean value indicating whether we complete the matching.
     */
    bool _step() {
        if (!_update_highway_list()) {
            return true;
        }
        int best_lane = _choose_best_highway();
        cost += (*highway_list)[best_lane].cost;

#ifdef DEBUG
        // update matched strings
        int distance = (*highway_list)[best_lane].starting_point + (*highway_list)[best_lane].first_hurdle -
                       (current_column + linear_leap_forward_column(current_lane, best_lane));
        _update_match(best_lane, current_lane, distance);
#endif
        // Update CIGAR
        _update_CIGAR(best_lane, current_lane, distance);
        CIGAR[CIGAR_index] = '\0';

        // Update position
        current_lane = best_lane;
        current_column = (*highway_list)[best_lane].starting_point + (*highway_list)[best_lane].first_hurdle;
#ifdef DEBUG
        printf("current position: %d, %d\n", current_lane, current_column);
#endif
        // Check if we reach the destination
        if (current_column >= (*highway_list)[current_lane].destination) {
            return true;
        }
        return false;
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
            lanes_without_short_hurdles[lane + k] = mask.flip_short_hurdles();
        }
    }

public:
    int_128bit& operator[](int lane){
        return lanes[lane + k];
    }

    hurdle_matrix(const char* read, const char* ref, int error) {
        m = std::min((size_t) MAX_LENGTH, strlen(read));
        n = std::min((size_t) MAX_LENGTH, strlen(ref));

        // assign to class parameters
        strncpy(A, read, m);
        strncpy(B, ref, n);
        k = std::max(error, abs(m-n) + 2);
        _convert_read();
        highway_list = new highways(k, m, n);
        lanes = new int_128bit[2 * k + 1];
        lanes_without_short_hurdles = new int_128bit[2 * k + 1];
        destination_lane = n - m;
        _construct_hurdles();

        // define starting position at (0, 0)
        current_lane = 0;
        current_column = 0;
        cost = 0;

#ifdef DEBUG
        strncpy(A_orig, read, m);
        strncpy(B_orig, ref, n);
        A_index = 0, B_index = 0, A_match_index = 0, B_match_index = 0;
#endif
        // initialize CIGAR string
        //CIGAR = new char[MAX_LENGTH * 2];
        CIGAR_index = 0;
    }

    void run() {
        bool flag = false;
        while (!flag) {
            flag = _step();
        }
        // Check if we reach the final destination
        int destination_column = (*highway_list)[destination_lane].destination;
        if (current_lane != destination_lane || current_column < destination_column) {
            printf("not reaching destination! current %d, %d -> %d; destination %d, %d\n",
                   current_lane, current_column, (*highway_list)[current_lane].destination, destination_lane, destination_column);
            int leap_cost = linear_leap_lane_penalty(current_lane, destination_lane);
            int hurdle_cost = std::max(0, destination_column - current_column -
                                       linear_leap_forward_column(current_lane, destination_lane));
            cost += leap_cost + hurdle_cost;
#ifdef DEBUG
            // update matched strings
            _update_match(destination_lane, current_lane, hurdle_cost);
#endif
            // update CIGAR string
            _update_CIGAR(destination_lane, current_lane, hurdle_cost);
        }
        printf("total cost: %d", cost);
    }

    void print() {
        for (int i = -k; i <= k; i++) {
            (*this)[i].print();
        }
        for (int i = -k; i <= k; i++) {
            this->lanes_without_short_hurdles[i + k].print();
        }
    }

    std::string get_CIGAR() {
        return CIGAR;
    }

    ~hurdle_matrix() {
        delete highway_list;
        delete[] lanes;
        delete[] lanes_without_short_hurdles;
    }
};

#endif //GASMA_UTILS_H
