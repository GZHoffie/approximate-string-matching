//
// Created by Zhenhao on 11/12/2021.
//

#ifndef GASMA_HURDLE_MATRIX_FLIPPED_H
#define GASMA_HURDLE_MATRIX_FLIPPED_H

#include "utils.h"
#include "hurdle_matrix.h"

class hurdle_matrix_flipped : public hurdle_matrix {
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
        int cost;

        // The starting column of the long highway
        int starting_point_long;

        // Length of the long highway
        int length_long;

        // Cost to reach the long highway
        int cost_long;

        // Destination column
        int destination;
    };

    // list containing the closest highway of each lane
    class highways {
    private:
        // maximum difference allowed
        int k;
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
        explicit highways(int error, int m, int n) {
            k = error;
            info = new highway_info[2 * k + 1];
            best_highway_lane = 0;
            for (int lane = -k; lane <= k; lane++) {
                info[lane + k].starting_point = -1;
                info[lane + k].length = 0;
                info[lane + k].cost = MAX_LENGTH;

                info[lane + k].starting_point_long = -1;
                info[lane + k].length_long = 0;
                info[lane + k].cost_long = MAX_LENGTH;

                info[lane + k].destination = _calculate_destination(m, n, lane);
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


protected:
    // rows in the hurdle matrix
    int_128bit* lanes_long;

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

            if ((*highway_list)[lane].starting_point_long < start_col) {
                // get closest highway in the lane
                int_128bit l = (lanes_long[lane + k]).shift_left(start_col);

                // update highway in lane
                first_zero = l.first_zero();
                int next_hurdle = (l.shift_left(first_zero)).first_one();
                (*highway_list)[lane].starting_point_long = start_col + first_zero;
                (*highway_list)[lane].length_long = next_hurdle;

                // Fix length if reaches destination
                if (start_col + first_zero + next_hurdle > (*highway_list)[lane].destination) {
                    (*highway_list)[lane].length_long = std::max(0, (*highway_list)[lane].destination -\
                                                   (start_col + first_zero));
                    reaching_destination = true;
                }
            }
            // calculate cost to reach the highway
            // FIXME: use function pointer for more complicated penalties.
            int leap_cost = linear_leap_lane_penalty(current_lane, lane);
            int hurdle_cost = (*highway_list)[lane].starting_point - start_col;
            (*highway_list)[lane].cost = leap_cost + hurdle_cost;

            hurdle_cost = (*highway_list)[lane].starting_point_long - start_col;
            (*highway_list)[lane].cost_long = leap_cost + hurdle_cost;


        }
        for (int lane = -k; lane <= k; lane++) {
            // get the best-looking highway
            int heuristic = std::min((*highway_list)[lane].cost_long - (*highway_list)[lane].length_long,
                                     (*highway_list)[lane].cost - (*highway_list)[lane].length);
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
    int _choose_best_highway() override {
        int best_lane = highway_list->best_highway_lane;
        int starting_point = (*highway_list)[best_lane].starting_point;
        int smallest_cost = (*highway_list)[best_lane].cost;
        int smallest_cost_length = 0;
        int smallest_cost_lane = best_lane;
        for (int lane = -k; lane <= k; lane++) {
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
        return smallest_cost_lane;
    }

    /**
     * Perform one step in the greedy algorithm.
     * @return a boolean value indicating whether we complete the matching.
     */
    bool _step() override{
        if (!_update_highway_list()) {
            return true;
        }
        int best_lane = _choose_best_highway();
        cost += (*highway_list)[best_lane].cost;

#ifdef DEBUG
        // update matched strings
        int distance = (*highway_list)[best_lane].starting_point + (*highway_list)[best_lane].length -
                       (current_column + linear_leap_forward_column(current_lane, best_lane));
        _update_match(best_lane, current_lane, distance);
#endif
        // Update CIGAR
        _update_CIGAR(best_lane, current_lane, distance);
        CIGAR[CIGAR_index] = '\0';

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

public:
    hurdle_matrix_flipped(const char *read, const char *ref,
                          int error) : hurdle_matrix(read, ref, error){
        m = std::min((size_t) MAX_LENGTH, strlen(read));
        n = std::min((size_t) MAX_LENGTH, strlen(ref));

        // assign to class parameters
        strncpy(A, read, m);
        strncpy(B, ref, n);
        k = std::max(error, abs(m-n) + 2);
        _convert_read();
        highway_list = new highways(k, m, n);
        lanes = new int_128bit[2 * k + 1];
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

        // construct long highways matrix
        lanes_long = new int_128bit[2 * k + 1];
        for (int lane = -k; lane <= k; lane++) {
            lanes_long[lane + k] = lanes[lane + k].flip_short_matches().flip_short_hurdles();
        }
    }

    void run() override {
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


    void print() override{
        for (int i = -k; i <= k; i++) {
            (*this)[i].print();
        }
        for (int i = -k; i <= k; i++) {
            lanes_long[i + k].print();
        }
    }

    ~hurdle_matrix_flipped() {
        delete[] lanes_long;
    }



};

#endif //GASMA_HURDLE_MATRIX_FLIPPED_H
