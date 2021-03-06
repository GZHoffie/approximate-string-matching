//
// Created by Zhenhao on 11/12/2021.
//

#ifndef GASMA_HURDLE_MATRIX_H
#define GASMA_HURDLE_MATRIX_H

#define MAX_K 50  // The maximum probable value for k

#include "utils.h"
#include <cstdlib>
#include <limits>
#include <math.h>

/**
 * The main class for the greedy string matching algorithm.
 * @tparam T either `int_128bit` (SSE) or `int_256bit` (AVX2), representing the type to store
 * the hurdle matrix in bits.
 */
template <typename T>
class hurdle_matrix {
private:
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
        int switch_cost;
        int hurdle_cost;

        // Number of actions to reach this highway
        int num_switches;
        int num_hurdles;

        // Destination column
        int destination;
    };

    // list containing the closest highway of each lane
    class highways {
    private:
        // maximum difference allowed
        int k;

        // upper and lower bound for the lane index
        int lower_bound, upper_bound;

        // highway information
        highway_info* info;

        static int _calculate_destination(int _m, int _n, int lane) {
            if (_m >= _n) {
                if (lane > 0) return _n - lane;
                else if (lane >= _n - _m) return _n;
                else return _m + lane;
            } else {
                if (lane < 0) return _m + lane;
                else if (lane <= _n - _m) return _m;
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
        explicit highways(int error, int m, int n, int lower_bound_, int upper_bound_) {
            k = error;
            lower_bound = lower_bound_;
            upper_bound = upper_bound_;
            info = new highway_info[2 * MAX_K + 1];
            best_highway_lane = 0;
            for (int lane = - MAX_K; lane <= MAX_K; lane++) {
                info[lane + MAX_K].starting_point = -1;
                info[lane + MAX_K].switch_cost = MAX_LENGTH;
                info[lane + MAX_K].hurdle_cost = MAX_LENGTH;
                info[lane + MAX_K].num_switches = MAX_LENGTH;
                info[lane + MAX_K].num_hurdles = MAX_LENGTH;
                info[lane + MAX_K].destination = _calculate_destination(m, n, lane);
            }
        }

        highway_info& operator[](int lane) {
            return info[lane + MAX_K];
        }

        void print() {
            for (int lane = lower_bound; lane <= upper_bound; lane++) {
                printf("lane: %d, starting point: %d, length: %d, cost: %d\n", lane,
                       (*this)[lane].starting_point, (*this)[lane].length, (*this)[lane].switch_cost + (*this)[lane].hurdle_cost);
            }
        }

        void reset(int error, int m_, int n_, int lower_bound_, int upper_bound_) {
            k = error;
            lower_bound = lower_bound_;
            upper_bound = upper_bound_;
            best_highway_lane = 0;
            for (int lane = lower_bound; lane <= upper_bound; lane++) {
                info[lane + MAX_K].starting_point = -1;
                info[lane + MAX_K].switch_cost = MAX_LENGTH;
                info[lane + MAX_K].hurdle_cost = MAX_LENGTH;
                info[lane + MAX_K].num_switches = MAX_LENGTH;
                info[lane + MAX_K].num_hurdles = MAX_LENGTH;
                info[lane + MAX_K].destination = _calculate_destination(m_, n_, lane);
            }
        }

        ~highways(){
            delete[] info;
        }
    };
    highways* highway_list;


protected:
    // maximum difference allowed
    int k;

    // upper and lower bound for the lane index
    int lower_bound, upper_bound;

    // two strings for comparison
    char A[MAX_LENGTH] __aligned;
    char B[MAX_LENGTH] __aligned;

#ifdef DISPLAY
    // string storing the original two strings
    char A_orig[MAX_LENGTH] __aligned;
    char B_orig[MAX_LENGTH] __aligned;

    // strings storing the matched strings
    char A_match[MAX_LENGTH * 2] __aligned;
    char B_match[MAX_LENGTH * 2] __aligned;
    int A_index, B_index, A_match_index, B_match_index;
#endif

    // T objects storing the bit arrays converted from string A and B
    T* A_bit0_mask;
    T* A_bit1_mask;
    T* B_bit0_mask;
    T* B_bit1_mask;

    // rows in the hurdle matrix
    T* lanes;

    // original rows (without flipping hurdles)
    T* lanes_orig;

    // information about destination
    int destination_lane;

    // length of read and ref
    int m, n;

    // CIGAR string used in mapper
    std::string CIGAR;

    // current position
    int current_lane;
    int current_column;

    // total cost
    int cost;

    // type of alignment
    alignment_type_t alignment_type;

    // mismatch penalty
    int x;

    // gap opening penalty
    int o;

    // gap extension penalty
    int e;

    // boolean value indicating whether it is the first step
    bool is_first_step;

    // significance calculation
    double match_sig, mismatch_sig, indel_sig;

#ifdef DISPLAY
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
#ifdef DEBUG
        printf("%s\n%s\n", A_match, B_match);
#endif
    }
#endif

    /**
     * Update the CIGAR strings given the lanes we are leaping to and the distance we
     * are moving.
     * @param best_lane The lane we are leaping to.
     * @param curr_lane The lane we are currently at.
     * @param mismatches The number of mismatches we need to go through before entering highway.
     * @param matches The length of highway we are going through.
     */
    void _update_CIGAR(int best_lane, int curr_lane, int mismatches, int matches) {
        // update matched strings
        if (best_lane < curr_lane) {
            CIGAR += std::to_string(curr_lane - best_lane);
            CIGAR += 'I';
        } else if (best_lane > curr_lane) {
            CIGAR += std::to_string(best_lane - curr_lane);
            CIGAR += 'D';
        }
        if (mismatches + matches > 0) {
            CIGAR += std::to_string(mismatches + matches);
            CIGAR += 'M';
        }
    }




    /**
     * Convert A and B into int8 objects and store them in A_bit0_t, A_bit1_t,
     * B_bit0_t and B_bit1_t.
     */
    void _convert_read() {
        // array of int8 objects to store the converted bits
        uint8_t A_bit0_t[MAX_LENGTH / 8] __aligned;
        uint8_t A_bit1_t[MAX_LENGTH / 8] __aligned;
        uint8_t B_bit0_t[MAX_LENGTH / 8] __aligned;
        uint8_t B_bit1_t[MAX_LENGTH / 8] __aligned;

        // convert string A and B into bits and store in the int8 array
        sse3_convert2bit1(A, A_bit0_t, A_bit1_t);
        sse3_convert2bit1(B, B_bit0_t, B_bit1_t);

        // convert the int8 array into 128-bit integer
        A_bit0_mask = new T(A_bit0_t);
        A_bit1_mask = new T(A_bit1_t);
        B_bit0_mask = new T(B_bit0_t);
        B_bit1_mask = new T(B_bit1_t);
    }


    /**
     * Update highway_list for each lane and show the closest highway to the
     * current position. Calculate the cost and record the longest highway.
     *
     * @return a boolean value indicating whether there is still a valid highway.
     */
    bool _update_highway_list() {
        double largest_total_heuristic = - std::numeric_limits<double>::infinity();
        int largest_leap_heuristic = - std::numeric_limits<int>::infinity();
        int best_highway_lane = 0;
        int first_zero;
        bool reaching_destination = false; // check if we are reaching destination
        for (int lane = lower_bound; lane <= upper_bound; lane++) {
            int start_col = current_column + switch_forward_column(current_lane, lane);
            if ((*highway_list)[lane].starting_point < start_col) {
                (*highway_list)[lane].num_switches = abs(lane - current_lane);
                // get closest highway in the lane
                T l = (lanes[lane + MAX_K]).shift_left(start_col);

                // update highway in lane
                first_zero = l.first_zero();
                int next_hurdle = (l.shift_left(first_zero)).first_one();
                (*highway_list)[lane].starting_point = start_col + first_zero;
                (*highway_list)[lane].length = next_hurdle;

                // Fix length if reaches destination
                if (start_col + first_zero + next_hurdle > (*highway_list)[lane].destination) {
                    (*highway_list)[lane].length = std::max(0, (*highway_list)[lane].destination -\
                                                   (start_col + first_zero));
                    reaching_destination = true;
                }
            }
            // calculate cost to reach the highway
            // FIXME: use function pointer for more complicated penalties.
            int switch_cost = 0;
            if (alignment_type == GLOBAL || !is_first_step) {
                switch_cost = switch_lane_penalty(current_lane, lane, o, e);
            }
            (*highway_list)[lane].num_hurdles = lanes_orig[lane + MAX_K].pop_count_between(start_col,(*highway_list)[lane].starting_point + (*highway_list)[lane].length);
            int hurdle_cost = x * (*highway_list)[lane].num_hurdles;
            (*highway_list)[lane].switch_cost = switch_cost;
            (*highway_list)[lane].hurdle_cost = hurdle_cost;

        }
        double heuristic;
        int leap_heuristic;
        for (int lane = lower_bound; lane <= upper_bound; lane++) {
            // get the best-looking highway
            int current_cost = - (*highway_list)[lane].switch_cost - (*highway_list)[lane].hurdle_cost;
            double significance = match_sig * (*highway_list)[lane].length +
                                  mismatch_sig * (*highway_list)[lane].num_hurdles +
                                  indel_sig * (*highway_list)[lane].num_switches;
            heuristic = significance;
            leap_heuristic = - (*highway_list)[lane].switch_cost;

            if (reaching_destination) {
                int final_switch_cost = 0;
                if (alignment_type == GLOBAL) {
                    final_switch_cost = switch_lane_penalty(lane, destination_lane, o, e);
                }
                heuristic = current_cost - final_switch_cost - x * ((*highway_list)[lane].destination -
                        (*highway_list)[lane].starting_point - (*highway_list)[lane].length);
                leap_heuristic -= final_switch_cost;
            }

            //printf("lane %d: heuristic %d\n", lane, heuristic);
            if (heuristic > largest_total_heuristic || (
                    heuristic == largest_total_heuristic && leap_heuristic > largest_leap_heuristic
                    )) {
                largest_total_heuristic = heuristic;
                largest_leap_heuristic = leap_heuristic;
                best_highway_lane = lane;
            }
        }
        highway_list->best_highway_lane = best_highway_lane;
#ifdef DEBUG
        highway_list->print();
        printf("Best highway lane: %d\n", best_highway_lane);
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
        // information about the highway on the best lane
        int best_lane = highway_list->best_highway_lane;
        int starting_point = (*highway_list)[best_lane].starting_point;
        int best_lane_cost = (*highway_list)[best_lane].hurdle_cost + (*highway_list)[best_lane].switch_cost;

        // keep the intermediate highway of smallest cost
        int smallest_intermediate_cost = best_lane_cost;
        int smallest_total_cost = best_lane_cost;
        int best_intermediate_lane = best_lane;

        // check all the other lanes for better highway
        best_lane = best_intermediate_lane;
        int intermediate_cost, total_cost, ending_point;
        for (int lane = lower_bound; lane <= upper_bound; lane++) {
            if (lane != best_lane) {
                if ((*highway_list)[lane].starting_point + switch_forward_column(lane, best_lane) > starting_point) {
                    continue;
                }
                ending_point = (*highway_list)[lane].starting_point + (*highway_list)[lane].length;
                intermediate_cost = (*highway_list)[lane].switch_cost + lanes_orig[lane + MAX_K].pop_count_between(current_column + switch_forward_column(current_lane, lane), ending_point);
                total_cost = intermediate_cost + switch_lane_penalty(lane, best_lane, o, e)
                             + std::max(0, x * lanes_orig[best_lane + MAX_K].pop_count_between(switch_forward_column(lane, best_lane) + ending_point, starting_point));
                if (total_cost <= smallest_total_cost) {
                    if (intermediate_cost <= smallest_intermediate_cost) {
                        smallest_total_cost = total_cost;
                        smallest_intermediate_cost = intermediate_cost;
                        best_intermediate_lane = lane;
                    }
                }
            }
        }
        return best_intermediate_lane;
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
        cost += (*highway_list)[best_lane].switch_cost + (*highway_list)[best_lane].hurdle_cost;

        // update matched strings
        int distance = (*highway_list)[best_lane].starting_point + (*highway_list)[best_lane].length -
                       (current_column + switch_forward_column(current_lane, best_lane));
#ifdef DISPLAY
        _update_match(best_lane, current_lane, distance);
#endif
        // Update CIGAR
        _update_CIGAR(best_lane, current_lane, distance - (*highway_list)[best_lane].length, (*highway_list)[best_lane].length);

        // Update position
        current_lane = best_lane;
        current_column = (*highway_list)[best_lane].starting_point + (*highway_list)[best_lane].length;
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
        T mask_bit0, mask_bit1;
        for (int lane = lower_bound; lane <= upper_bound; lane++) {
            if (lane < 0) {
                mask_bit0 = (A_bit0_mask->shift_left(-lane))._xor(*B_bit0_mask);
                mask_bit1 = (A_bit1_mask->shift_left(-lane))._xor(*B_bit1_mask);
            } else {
                mask_bit0 = (B_bit0_mask->shift_left(lane))._xor(*A_bit0_mask);
                mask_bit1 = (B_bit1_mask->shift_left(lane))._xor(*A_bit1_mask);
            }
            auto mask = mask_bit0._or(mask_bit1);
            lanes_orig[lane + MAX_K] = mask;
            lanes[lane + MAX_K] = mask.flip_short_hurdles(1);//.flip_short_matches(1);
        }
    }

public:
    T& operator[](int lane){
        return lanes[lane + MAX_K];
    }

    /**
     * Default Constructor of the hurdle matrix
     * @param read the read string
     * @param ref the reference string
     * @param error the maximum consecutive insert/delete allowed
     * @param _alignment_type the type of alignment, either GLOBAL, SEMI_GLOBAL or LOCAL (currently not supported).
     *                       Default: GLOBAL.
     * @param _x penalty for mismatch. Default: 1.
     * @param _o gap opening penalty. Default: 0.
     * @param _e gap extension penalty. Default: 1.
     */
    hurdle_matrix(
            const char* read,
            const char* ref,
            int error,
            alignment_type_t _alignment_type = GLOBAL,
            int _x = 1,
            int _o = 1,
            int _e = 1,
            double match_prob = 0.95,
            double mismatch_prob = 0.02,
            double indel_prob = 0.03
            ) {
        // TODO: replace MAX_LENGTH with either 128 or 256, depending on the whether
        // `T` is `int_128bit` or `int_256bit`.
        m = std::min(MAX_LENGTH, static_cast<int>(strlen(read)));
        n = std::min(MAX_LENGTH, static_cast<int>(strlen(ref)));

        // Set alignment parameters
        alignment_type = _alignment_type;
        x = _x;
        o = _o;
        e = _e;

        // assign to class parameters
        strncpy(A, read, m);
        strncpy(B, ref, n);
        k = error;

#ifdef CORRECTION
        if (m <= n) {
            lower_bound = -k;
            upper_bound = n - m + k;
        } else {
            lower_bound = n - m - k;
            upper_bound = k;
        }
#else
        lower_bound = -k;
        upper_bound = k;
#endif

        _convert_read();
        highway_list = new highways(MAX_K, m, n, lower_bound, upper_bound);
        lanes = new T[2 * MAX_K + 1];
        lanes_orig = new T[2 * MAX_K + 1];
        destination_lane = n - m;
        is_first_step = true;
        _construct_hurdles();

        // define starting position at (0, 0)
        current_lane = 0;
        current_column = 0;
        cost = 0;

#ifdef DISPLAY
        strncpy(A_orig, read, m);
        strncpy(B_orig, ref, n);
        A_index = 0, B_index = 0, A_match_index = 0, B_match_index = 0;
#endif
        // initialize CIGAR string
        CIGAR.reserve(MAX_LENGTH * sizeof(char));

        // calculate significance for match/mismatch/indel
        match_sig = log(match_prob / 0.25);
        mismatch_sig = log(mismatch_prob / 0.25);
        indel_sig = log(indel_prob / 2 / 0.25);
    }



    /**
     * Constructor of the class that sets alignment type and penalty scheme only.
     * Default: use the edit distance penalty scheme.
     * @param _alignment_type the type of alignment, either GLOBAL, SEMI_GLOBAL or LOCAL (currently not supported).
     *                       Default: GLOBAL.
     * @param _x penalty for mismatch. Default: 1.
     * @param _o gap opening penalty. Default: 0.
     * @param _e gap extension penalty. Default: 1.
     */
    explicit hurdle_matrix(
            alignment_type_t _alignment_type = GLOBAL,
            int _x = 1,
            int _o = 1,
            int _e = 1,
            double match_prob = 0.80,
            double mismatch_prob = 0.20 / 3,
            double indel_prob = 0.40 / 3
                    ) : hurdle_matrix("", "", 1, _alignment_type, _x, _o, _e,
                                      match_prob, mismatch_prob, indel_prob)
                    {}


    /**
     * Run the greedy algorithm for several steps.
     */
    void run() {
        bool flag = false;
        while (!flag) {
            flag = _step();
            is_first_step = false;
        }
        // Check if we reach the final destination
        int destination_column = (*highway_list)[destination_lane].destination;
        if (current_lane != destination_lane || current_column < destination_column) {
            int switch_cost = 0;
            if (alignment_type == GLOBAL) {
                switch_cost = switch_lane_penalty(current_lane, destination_lane, o, e);
            }
            int distance = lanes_orig[destination_lane + MAX_K].pop_count_between(current_column + switch_forward_column(current_lane, destination_lane), destination_column);
            int hurdle_cost = std::max(0, x * distance);
            cost += switch_cost + hurdle_cost;
#ifdef DISPLAY
            // update matched strings
            _update_match(destination_lane, current_lane, hurdle_cost);
#endif
            // update CIGAR string
            _update_CIGAR(destination_lane, current_lane, distance, 0);
        }
#ifdef DISPLAY
        A_match[A_match_index] = '\0';
        B_match[B_match_index] = '\0';
        printf("%s\n%s\n", A_match, B_match);
#endif
        //printf("total cost: %d", cost);
    }

    /**
     * Print out the hurdle matrix in bit form.
     */
    void print() {
        for (int i = lower_bound; i <= upper_bound; i++) {
            printf("lane %d:", i);
            (*this)[i].print();
        }
    }

    /**
     * Get the CIGAR string. Must be called after run().
     * @return CIGAR string.
     */
    std::string get_CIGAR() {
        return CIGAR;
    }

    /**
     * Reset the object to get ready for the next alignment.
     * @param read the read string.
     * @param read_len length of the read string.
     * @param ref the reference string.
     * @param ref_len length of the reference string.
     * @param error band width
     */
    void reset(const char* read, const int read_len, const char* ref, const int ref_len, int error) {
        m = std::min(MAX_LENGTH, read_len);
        n = std::min(MAX_LENGTH, ref_len);

        // assign to class parameters
        strncpy(A, read, m);
        strncpy(B, ref, n);
        k = error;

#ifdef CORRECTION
        if (m <= n) {
            lower_bound = -k;
            upper_bound = n - m + k;
        } else {
            lower_bound = n - m - k;
            upper_bound = k;
        }
#else
        lower_bound = -k;
        upper_bound = k;
#endif

        _convert_read();
        highway_list->reset(k, m, n, lower_bound, upper_bound);
        destination_lane = n - m;
        is_first_step = true;
        _construct_hurdles();

        // define starting position at (0, 0)
        current_lane = 0;
        current_column = 0;
        cost = 0;

#ifdef DISPLAY
        strncpy(A_orig, read, m);
        strncpy(B_orig, ref, n);
        A_index = 0, B_index = 0, A_match_index = 0, B_match_index = 0;
#endif
        // initialize CIGAR string
        CIGAR.clear();
    }

    void reset(const char* read, const char* ref, int error) {
        int read_len = static_cast<int>(strlen(read));
        int ref_len = static_cast<int>(strlen(ref));
        reset(read, read_len, ref, ref_len, error);
    }

    /**
     * Return the penalty (non-negative) of the alignment.
     * @return total penalty
     */
    int get_cost() const {
        return cost;
    }

    ~hurdle_matrix() {
        delete highway_list;
        delete[] lanes;
    }
};




#endif //GASMA_HURDLE_MATRIX_H
