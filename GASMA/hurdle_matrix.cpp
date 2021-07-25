//
// Created by 顾振昊 on 2021/5/22.
//

#include "hurdle_matrix.h"
#include <x86intrin.h>

hurdle_matrix::hurdle_matrix(int k, bool remove_single_zeros) {
    this->k = k;
    this->remove_single_zeros = remove_single_zeros;
    this->mat = new __m256i[2 * this->k + 1];
    this->initialized = true;
}

hurdle_matrix::~hurdle_matrix() {
    delete[] mat;
}

void hurdle_matrix::load_reads(const std::string &read, const std::string &ref) {
    for (int i = -k; i <= k; i++) {
        this->mat[i + k] = this->hurdle_info(read, ref, i);
    }
}

__m256i hurdle_matrix::hurdle_info(const std::string &read, const std::string &ref, int shift) {

    return 0;
}
