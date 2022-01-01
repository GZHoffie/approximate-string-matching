//
// Created by Zhenhao on 2/1/2022.
//

#ifndef LEAP_SIMD_SIMD_ED_H
#define LEAP_SIMD_SIMD_ED_H


#include <iostream>
#include <x86intrin.h>
#include <cstring>
#include <string>
#include <cstdlib>
#include "print.h"
#include "shift.h"
#include "bit_convert.h"

#ifndef _MAX_LENGTH_
#define _MAX_LENGTH_ 256
#endif

#ifndef _AFFINE_DEF_POS_
#define _AFFINE_DEF_POS_ 10000000
#endif

using namespace std;

#ifndef __ED_INFO_H_
#define __ED_INFO_H_

enum ED_TYPE {MISMATCH, A_INS, B_INS};

struct ED_INFO {
    ED_TYPE type;
    int id_length;
};

#endif


enum ED_modes {ED_LOCAL, ED_GLOBAL, ED_SEMI_FREE_BEGIN, ED_SEMI_FREE_END};
enum OP_modes {SSE, AVX};

class SIMD_ED {
public:
    SIMD_ED();
    ~SIMD_ED();

    void init_levenshtein(int ED_threshold, ED_modes mode = ED_LOCAL, bool SHD_enable = true);
    void init_affine(int gap_threshold, int AF_threshold, ED_modes mode, int ms_penalty, int gap_open_penalty, int gap_ext_penalty, bool SHD_enable = false, int SHD_threshold = 10);
    int count_ID_length_avx(int lane_idx, int start_pos);

    void convert_reads(char *read, char *ref, int length, uint8_t *A0, uint8_t *A1, uint8_t *B0, uint8_t *B1);

    void load_reads(char *read, char *ref, int length);
    void load_reads(uint8_t *A0, uint8_t *A1, uint8_t *B0, uint8_t *B1, int length);
    void load_reads(__m256i A0, __m256i A1, __m256i B0, __m256i B1, int length);

    void load_ref(__m256i B0, __m256i B1);
    void load_read(__m256i A0, __m256i A1, int length);

    void calculate_masks();

    void reset();
    void run();
    bool check_pass();
    void backtrack();
    int get_ED();
    string get_CIGAR();
private:
    void reset_levenshtein();
    void run_levenshtein();
    void backtrack_levenshtein();

    void run_affine();
    void reset_affine();
    void backtrack_affine();

    // variables
    int ED_t;
    __m256i *hamming_masks;
    ED_modes mode;
    bool SHD_enable;
    //__m128i shifted_mask;

    // Affine mode
    int SHD_threshold;
    bool affine_mode;
    int gap_threshold;
    int af_threshold;
    int ms_penalty;
    int gap_open_penalty;
    int gap_ext_penalty;
    int **I_pos;
    int **D_pos;
    int ED_count;

    // information of each lane
    int *cur_ED;
    int **start;
    int **end;

    // backtracking data
    bool ED_pass;
    int final_lane_idx;
    int final_ED;
    ED_INFO *ED_info;
    // data only for global and semi_global
    int converge_ED;
    int converge_final_lane;

    // serrogates
    int mid_lane;
    int total_lanes;

    // buffers
    int buffer_length;

    char A[_MAX_LENGTH_] __aligned;
    char B[_MAX_LENGTH_] __aligned;

    uint8_t A_bit0_t[_MAX_LENGTH_ / 4] __aligned;
    uint8_t A_bit1_t[_MAX_LENGTH_ / 4] __aligned;
    uint8_t B_bit0_t[_MAX_LENGTH_ / 4] __aligned;
    uint8_t B_bit1_t[_MAX_LENGTH_ / 4] __aligned;

};


#endif //LEAP_SIMD_SIMD_ED_H
