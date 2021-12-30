//
// Created by Zhenhao on 28/12/2021.
//

#ifndef GASMA_BENCHMARK_UTILS_H
#define GASMA_BENCHMARK_UTILS_H

#include <iostream>
#include <string>
#include <sys/times.h>
#include <unistd.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>

#include "LEAP/SIMD_ED.h"
#include "parasail/parasail.h"

#include "../hurdle_matrix.h"


/**
 * Class for benchmarking. Algorithms for benchmarking include
 * 1. Greedy
 * 2. Needleman-Wunsch (https://github.com/jeffdaily/parasail)
 * 3. LEAP (https://github.com/CMU-SAFARI/LEAP)
 */
class benchmark {
private:
    // LEAP objects
    SIMD_ED* ed_obj;

    // SIMD Needleman-Wunsch objects
    parasail_matrix_t* penalty_matrix;

    // Greedy algorithm objects
    hurdle_matrix* matrix;

    // whether we use SIMD acceleration for NW and LEAP
    bool use_SIMD;

    // mismatch penalty (non-negative)
    int x;

    // gap opening penalty (non-negative)
    int o;

    // gap extension penalty (non-negative)
    int e;

    // band width
    int k;

    // maximum test number
    int max_tests;

    // array of strings storing the strings for alignment
    std::string* read;
    std::string* ref;

    // array storing the optimal results
    int* answers;

    // tms objects to record the time usage
    tms start_time;
    tms end_time;
    tms nw_time, LEAP_time, greedy_time;

    // check correctness of the algorithm
    int total_tests;
    int nw_correct, LEAP_correct, greedy_correct;

    /**
     * Run banded Needleman-Wunsch algorithm on two strings s1 and s2.
     * @param s1 the read string.
     * @param s1Len length of the read string.
     * @param s2 the reference string.
     * @param s2Len length of the reference string.
     * @return the penalty of alignment (non-negative).
     */
    int _run_nw_banded(
            const char * s1,
            const int s1Len,
            const char * s2,
            const int s2Len) {
        parasail_result_t* result;
        // TODO: check what the band width is affecting
        times(&start_time);
        result = parasail_nw_banded(s1, s1Len, s2, s2Len, -o, -e, k, penalty_matrix);
        times(&end_time);
        nw_time.tms_stime += end_time.tms_stime - start_time.tms_stime;
        nw_time.tms_utime += end_time.tms_utime - start_time.tms_utime;
        return -result->score;
    }

    /**
     * Run Needleman-Wunsch algorithm with SSE acceleration on two strings s1 and s2.
     * @return the penalty of alignment (non-negative).
     */
    int _run_nw_sse(
            const char * s1,
            const int s1Len,
            const char * s2,
            const int s2Len) {
        parasail_result_t* result;
        times(&start_time);
        result =  parasail_nw_striped_sse41_128_16(s1, s1Len, s2, s2Len, -o, -e, penalty_matrix);
        times(&end_time);
        nw_time.tms_stime += end_time.tms_stime - start_time.tms_stime;
        nw_time.tms_utime += end_time.tms_utime - start_time.tms_utime;
        return -result->score;
    }

    /**
     * Run LEAP algorithm with SIMD acceleration on two strings s1 and s2.
     * @return the penalty of alignment (non-negative).
     */
    int _run_LEAP(
            char * s1,
            int s1Len,
            char * s2,
            int s2Len
            ) {
        int length = std::max(s1Len, s2Len);
        times(&start_time);
        ed_obj->load_reads(s1, s2, length);
        ed_obj->calculate_masks();
        ed_obj->reset();
        ed_obj->run();
        // TODO: check the efficiency of this backtrack step.
        if (ed_obj->check_pass()) {
            ed_obj->backtrack();
        }
        times(&end_time);
        LEAP_time.tms_stime += end_time.tms_stime - start_time.tms_stime;
        LEAP_time.tms_utime += end_time.tms_utime - start_time.tms_utime;
        return ed_obj->get_ED();
    }

    /**
     * Run greedy algorithm with SIMD acceleration on two strings s1 and s2.
     * @return the penalty of alignment (non-negative).
     */
    int _run_greedy(
            const char * s1,
            const int s1Len,
            const char * s2,
            const int s2Len
    ) {
        times(&start_time);
        matrix->reset(s1, s1Len, s2, s2Len, k);
        matrix->run();
        times(&end_time);
        greedy_time.tms_stime += end_time.tms_stime - start_time.tms_stime;
        greedy_time.tms_utime += end_time.tms_utime - start_time.tms_utime;
        return matrix->get_cost();
    }

    /**
     * Run benchmark on the three algorithms and record whether the results are correct.
     * @param correct_answer the correct cost of this alignment (non-negative).
     */
    void _run_benchmark(
            char * s1,
            int s1Len,
            char * s2,
            int s2Len,
            int correct_answer = -1
            ) {
        total_tests += 1;
        int nw_result, LEAP_result, greedy_result;
        if (use_SIMD) {
            nw_result = _run_nw_sse(s1, s1Len, s2, s2Len);
        } else {
            nw_result = _run_nw_banded(s1, s1Len, s2, s2Len);
        }
        LEAP_result = _run_LEAP(s1, s1Len, s2, s2Len);
        greedy_result = _run_greedy(s1, s1Len, s2, s2Len);

        // check correctness
        if (correct_answer == -1) {
            // set correct answer as the banded NW algorithm
            correct_answer = nw_result;
        }
        nw_correct += (nw_result == correct_answer);
        LEAP_correct += (LEAP_result == correct_answer);
        greedy_correct += (greedy_result == correct_answer);
    }


public:
    benchmark(
            int _x,
            int _o,
            int _e,
            int _k,
            int max_test_num,
            bool _use_SIMD = true)
            {
        // set scoring scheme
        x = _x;
        o = _o;
        e = _e;
        k = _k;

        // maximum alignment tests number
        max_tests = max_test_num;
        read = new std::string[max_tests];
        ref = new std::string[max_tests];
        answers = new int[max_tests];
        std::fill_n(answers, max_tests, -1);

        // initialize the objects used for benchmarking
        use_SIMD = _use_SIMD;
        ed_obj = new SIMD_ED;
        matrix = new hurdle_matrix(GLOBAL, x, o, e);
        penalty_matrix = parasail_matrix_create("ACGT", 0, -x);
        if (o == e && o == 1) {
            // Levenshtein distance
            // FIXME: here it seems that the edit distance is calculated
            ed_obj->init_levenshtein(k, ED_GLOBAL, use_SIMD);
        } else {
            // Affine gap penalty
            ed_obj->init_affine(k, k * 3, ED_GLOBAL, x, o, e, use_SIMD);
        }

        // Initialize time objects
        nw_time.tms_stime = 0;
        nw_time.tms_utime = 0;
        nw_time.tms_cstime = 0;
        nw_time.tms_cutime = 0;

        LEAP_time.tms_stime = 0;
        LEAP_time.tms_utime = 0;
        LEAP_time.tms_cstime = 0;
        LEAP_time.tms_cutime = 0;

        greedy_time.tms_stime = 0;
        greedy_time.tms_utime = 0;
        greedy_time.tms_cstime = 0;
        greedy_time.tms_cutime = 0;

        // initialize correctness record
        total_tests = 0;
        nw_correct = LEAP_correct = greedy_correct = 0;
    }

    /**
     * Read the file that contains read strings and reference strings,
     * each on a separate line.
     * @param string_dir the directory to the file
     * @param skip_first_char true if the first character on each line is to be skipped.
     */
    void read_string_file(
            const char * string_dir,
            bool skip_first_char = true) {
        std::ifstream string_file;
        string_file.open(string_dir);
        std::string read_temp, ref_temp;
        int i;
        for (i = 0; i < max_tests; i++) {
            if (skip_first_char) {
                string_file.ignore();
            }
            if (!getline(string_file, read_temp)) {
                break;
            }
            if (skip_first_char) {
                string_file.ignore();
            }
            getline(string_file, ref_temp);
            read[i].assign(read_temp);
            ref[i].assign(ref_temp);
        }
        max_tests = i;
    }

    /**
     * Read the file containing the corresponding optimal penalty scores.
     * @param answer_dir the directory to the file
     */
    void read_answer_file(const char * answer_dir) {
        std::ifstream answer_file;
        answer_file.open(answer_dir);
        std::string optimal_res_temp;
        for (int i = 0; i < max_tests; i++) {
            if (!getline(answer_file, optimal_res_temp)) {
                break;
            }
            answers[i] = atoi(optimal_res_temp.c_str());
        }
    }

    /**
     * Run the benchmark. Must be run after read_string_file().
     */
    void run() {
        for (int i = 0; i < max_tests; i++) {
            auto s1 = (char*) read[i].c_str();
            auto s2 = (char*) ref[i].c_str();
            auto s1Len = (int) read[i].length();
            auto s2Len = (int) ref[i].length();
            _run_benchmark(s1, s1Len, s2, s2Len, answers[i]);
            if (i % 10000 == 0 && i != 0) {
                printf("...processed %d reads.\n", i);
            }
        }
        printf("...complete.\n");
    }

    /**
     * Print out the benchmark results. Must be run after run().
     */
    void print() {
        printf("===================== Benchmark Results =====================\n");
        printf("Total number of alignments: %d\n[Time]\n", total_tests);
        printf("=> Needleman-Wunsch time: %f s\n", (double) nw_time.tms_utime / sysconf(_SC_CLK_TCK));
        printf("=> LEAP time: %f s\n", (double) LEAP_time.tms_utime / sysconf(_SC_CLK_TCK));
        printf("=> Greedy time: %f s\n", (double) greedy_time.tms_utime / sysconf(_SC_CLK_TCK));
        printf("[Accuracy]\n");
        printf("=> Needleman-Wunsch correct rate: %f\n", (double) nw_correct / total_tests);
        printf("=> LEAP correct rate: %f\n", (double) LEAP_correct / total_tests);
        printf("=> Greedy correct rate: %f\n", (double) greedy_correct / total_tests);
    }

    ~benchmark() {
        delete matrix;
        delete ed_obj;
        delete[] read;
        delete[] ref;
        delete[] answers;
    }
};


#endif //GASMA_BENCHMARK_UTILS_H
