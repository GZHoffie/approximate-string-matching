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


#include "parasail/parasail.h"
#include "../hurdle_matrix.h"
#include "LEAP_SIMD/LV_BAG.h"

#include "benchmark_coverage.h"

/**
 * Structure to store the results of alignment.
 * penalty: the alignment penalty.
 * CIGAR: CIGAR string storing the result of backtracking.
 */
class align_result_t {
public:
    int penalty;
    string CIGAR;

    align_result_t () {
        penalty = 0;
        CIGAR = "";
    }
};



/**
 * Class for benchmarking. Algorithms for benchmarking include
 * 1. Greedy
 * 2. Needleman-Wunsch (https://github.com/jeffdaily/parasail)
 * 3. LEAP (https://github.com/CMU-SAFARI/LEAP)
 */
class benchmark {
private:
    // LEAP objects
    LV* ed_obj;

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

    // align_result_t objects to store the alignment results
    align_result_t *nw_results, *LEAP_results, *greedy_results;

    // check correctness of the algorithm
    int total_tests;
    int nw_correct, LEAP_correct, greedy_correct;
    int greedy_coverage;

    /**
     * Run banded Needleman-Wunsch algorithm on two strings s1 and s2. Store the results
     * in nw_results.
     * @param s1 the read string.
     * @param s1Len length of the read string.
     * @param s2 the reference string.
     * @param s2Len length of the reference string.
     */
    void _run_nw(
            const char * s1,
            const int s1Len,
            const char * s2,
            const int s2Len) {
        parasail_result_t* result;
        parasail_cigar_t* cigar_result;
        // TODO: check what the band width is affecting
        times(&start_time);
        result = parasail_nw_trace(s1, s1Len, s2, s2Len, o, e, penalty_matrix);
        cigar_result = parasail_result_get_cigar(result, s1, s1Len, s2, s2Len, penalty_matrix);
        nw_results->CIGAR = parasail_cigar_decode(cigar_result);
        nw_results->penalty = - result->score;
        times(&end_time);

        nw_time.tms_stime += end_time.tms_stime - start_time.tms_stime;
        nw_time.tms_utime += end_time.tms_utime - start_time.tms_utime;

        parasail_result_free(result);
        parasail_cigar_free(cigar_result);
    }

    /**
     * Run Needleman-Wunsch algorithm with SSE acceleration on two strings s1 and s2. Store the results
     * in nw_results.
     */
    void _run_nw_sse(
            const char * s1,
            const int s1Len,
            const char * s2,
            const int s2Len) {
        parasail_result_t* result;
        parasail_cigar_t* cigar_result;

        times(&start_time);
        result =  parasail_nw_trace_striped_sse41_128_16(s1, s1Len, s2, s2Len, o, e, penalty_matrix);
        cigar_result = parasail_result_get_cigar(result, s1, s1Len, s2, s2Len, penalty_matrix);
        nw_results->CIGAR = parasail_cigar_decode(cigar_result);
        nw_results->penalty = - result->score;
        times(&end_time);

        nw_time.tms_stime += end_time.tms_stime - start_time.tms_stime;
        nw_time.tms_utime += end_time.tms_utime - start_time.tms_utime;

        parasail_result_free(result);
        parasail_cigar_free(cigar_result);
    }

    /**
     * Run LEAP algorithm with SIMD acceleration on two strings s1 and s2. Store the results
     * in LEAP_results.
     */
    void _run_LEAP(
            char * s1,
            int s1Len,
            char * s2,
            int s2Len
            ) {
        int length = std::max(s1Len, s2Len);

        times(&start_time);
        ed_obj->load_reads(s1, s2, length);
        //ed_obj->calculate_masks();
        ed_obj->reset();
        ed_obj->run();
        // TODO: check the efficiency of this backtrack step.
        if (ed_obj->check_pass()) {
            ed_obj->backtrack();
        }
        LEAP_results->penalty = ed_obj->get_ED();
        LEAP_results->CIGAR = ed_obj->get_CIGAR();
        times(&end_time);

        LEAP_time.tms_stime += end_time.tms_stime - start_time.tms_stime;
        LEAP_time.tms_utime += end_time.tms_utime - start_time.tms_utime;
    }

    /**
     * Run greedy algorithm on two strings s1 and s2.
     * Store the results in greedy_results.
     */
    void _run_greedy(
            const char * s1,
            const int s1Len,
            const char * s2,
            const int s2Len
    ) {
        times(&start_time);
        matrix->reset(s1, s1Len, s2, s2Len, k);
        matrix->run();
        greedy_results->penalty = matrix->get_cost();
        greedy_results->CIGAR = matrix->get_CIGAR();
        times(&end_time);

        greedy_time.tms_stime += end_time.tms_stime - start_time.tms_stime;
        greedy_time.tms_utime += end_time.tms_utime - start_time.tms_utime;
    }

    /**
     * Check if the LCM contained in the alignment 1 (as indicated in CIGAR1) covers that contained
     * in alignment 2 (as indicated in CIGAR2).
     * @param s1 the read string.
     * @param s2 the reference string.
     * @param CIGAR1 the CIGAR string given in the first alignment.
     * @param CIGAR2 the CIGAR string given in the second alignment.
     * @param threshold1 the threshold for consecutive matching in alignment 1 to be considered as LCM.
     * @param threshold2 the threshold for consecutive matching in alignment 2 to be considered as LCM.
     * @return a boolean value indicating whether alignment 1 covers alignment 2.
     */
    bool _check_coverage(
            char * s1,
            char * s2,
            string CIGAR1,
            string CIGAR2,
            int threshold1,
            int threshold2
            ) {
        std::string LCM1 = long_consecutive_matching_substring(s1, s2, CIGAR1, threshold1);
        std::string LCM2 = long_consecutive_matching_substring(s1, s2, CIGAR2, threshold2);
        return covers(LCM1, LCM2);
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
            int correct_answer = INT32_MIN
            ) {
        total_tests += 1;
        if (use_SIMD) {
            _run_nw_sse(s1, s1Len, s2, s2Len);
        } else {
            _run_nw(s1, s1Len, s2, s2Len);
        }
        _run_LEAP(s1, s1Len, s2, s2Len);
        _run_greedy(s1, s1Len, s2, s2Len);
        //printf("%s\n%s\n", s1, s2);
        //printf("%d, %d\n", nw_results->penalty, greedy_results->penalty);
        // check correctness
        if (correct_answer == INT32_MIN) {
            // set correct answer as the banded NW algorithm
            correct_answer = nw_results->penalty;
        }
        nw_correct += (nw_results->penalty == correct_answer);
        LEAP_correct += (LEAP_results->penalty == correct_answer);
        greedy_correct += (greedy_results->penalty == correct_answer);
        greedy_coverage += (int) (_check_coverage(s1, s2, greedy_results->CIGAR, nw_results->CIGAR, 1, 3));
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
        std::fill_n(answers, max_tests, INT32_MIN);

        // initialize the objects used for benchmarking
        use_SIMD = _use_SIMD;
        ed_obj = new LV;
        matrix = new hurdle_matrix(GLOBAL, x, o, e);
        penalty_matrix = parasail_matrix_create("ACGT", 0, -x);
        ed_obj->init(k, 200, ED_GLOBAL, x, o, e);

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

        // Initialize results
        nw_results = new align_result_t;
        LEAP_results = new align_result_t;
        greedy_results = new align_result_t;

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
        printf("=> Needleman-Wunsch | %.3f s\n", (double) nw_time.tms_utime / sysconf(_SC_CLK_TCK));
        printf("=> LEAP             | %.3f s\n", (double) LEAP_time.tms_utime / sysconf(_SC_CLK_TCK));
        printf("=> Greedy           | %.3f s\n", (double) greedy_time.tms_utime / sysconf(_SC_CLK_TCK));
        printf("[Accuracy] (percentage of alignments matching optimal penalty)\n");
        printf("=> Needleman-Wunsch | %.3f %%\n", (double) nw_correct / total_tests * 100);
        printf("=> LEAP             | %.3f %%\n", (double) LEAP_correct / total_tests * 100);
        printf("=> Greedy           | %.3f %%\n", (double) greedy_correct / total_tests * 100);
        printf("[Coverage] (percentage of alignments covering all long consecutive matches)\n");
        printf("=> Greedy           | %.3f %%\n", (double) greedy_coverage / total_tests * 100);
    }

    ~benchmark() {
        delete matrix;
        delete ed_obj;
        delete nw_results;
        delete LEAP_results;
        delete greedy_results;
        delete[] read;
        delete[] ref;
        delete[] answers;
    }
};


#endif //GASMA_BENCHMARK_UTILS_H
