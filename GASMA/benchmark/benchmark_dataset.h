//
// Created by Zhenhao on 3/1/2022.
//

#ifndef GASMA_BENCHMARK_DATASET_H
#define GASMA_BENCHMARK_DATASET_H

/*
 *                             The MIT License
 *
 * Wavefront Alignments Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignments Algorithms.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * PROJECT: Wavefront Alignments Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Sequence Generator for benchmarking pairwise algorithms
 */

#include <cstdio>
#include <cstdlib>
#include <getopt.h>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <ctime>

#ifndef ALPHABET_SIZE
#define ALPHABET_SIZE 4
#endif
char alphabet[] = {
        'A','C','G','T'
};


/**
 * Utility class to generate simulated dataset of read and ref strings
 * with certain amount of errors. This class is adopted from
 * https://github.com/smarco/WFA/blob/master/tools/generate_dataset.c
 */
class Dataset {
private:
    // number of reads to be output in the dataset
    int num_reads;

    // length of the read string
    int length;

    // error rate
    float error_rate;

    // whether the generated number of error is exactly equal to error_rate * length or <= that number.
    bool exact_error_rate;

    // whether we overwrite file if file already exist
    bool overwrite;

    // probability of mismatches
    float mismatch_rate;

    /**
     * Generate random number between min and max
     * @return a uint64 object that is randomly generated between min and max
     */
    uint64_t rand_iid(const uint64_t min,const uint64_t max) {
        int n_rand = rand(); // [0, RAND_MAX]
        const uint64_t range = max - min;
        const uint64_t rem = RAND_MAX % range;
        const uint64_t sample = RAND_MAX / range;
        // Consider the small interval within remainder of RAND_MAX
        if (n_rand < RAND_MAX - rem) {
            return min + n_rand/sample;
        } else {
            return rand_iid(min,max);
        }
    }
/*
 * Generate pattern
 */
    void generate_pattern(
            char* const pattern,
            const uint64_t length) {
        // Generate random characters
        uint64_t i;
        for (i=0;i<length;++i) {
            pattern[i] = alphabet[rand_iid(0,ALPHABET_SIZE)];
        }
        pattern[length] = '\0';
    }
/*
 * Generate candidate-text from pattern adding random errors
 */
    void generate_candidate_text_add_mismatch(
            char* const candidate_text,
            const uint64_t candidate_length) {
        // Generate random mismatch
        int position = rand_iid(0,candidate_length);
        char character = alphabet[rand_iid(0,ALPHABET_SIZE)];
        candidate_text[position] = character;
    }
    void generate_candidate_text_add_deletion(
            char* const candidate_text,
            uint64_t* const candidate_length) {
        // Generate random deletion
        int position = rand_iid(0,*candidate_length);
        int i;
        const uint64_t new_candidate_length = *candidate_length - 1;
        for (i=position;i<new_candidate_length;++i) {
            candidate_text[i] = candidate_text[i+1];
        }
        *candidate_length = new_candidate_length;
    }
    void generate_candidate_text_add_insertion(
            char* const candidate_text,
            uint64_t* const candidate_length) {
        // Generate random insertion
        int position = rand_iid(0,*candidate_length);
        int i;
        const uint64_t new_candidate_length = *candidate_length + 1;
        for (i=new_candidate_length-1;i>position;--i) {
            candidate_text[i] = candidate_text[i-1];
        }
        *candidate_length = new_candidate_length;
        // Insert random character
        candidate_text[position] = alphabet[rand_iid(0,ALPHABET_SIZE)];
    }
    char* generate_candidate_text(
            char* const pattern,
            const uint64_t pattern_length,
            const float error_degree) {
        // Compute nominal number of errors
        uint64_t num_errors;
        if (exact_error_rate)
            num_errors = ceil(pattern_length * error_degree);
        else
            num_errors = rand_iid(0, ceil(pattern_length * error_degree));
        // Allocate & init-by-copy candidate text
        char* const candidate_text = new char[pattern_length + num_errors];
        uint64_t candidate_length = pattern_length;
        memcpy(candidate_text,pattern,pattern_length);
        // Generate random errors
        int i;
        for (i=0;i<num_errors;++i) {
            float random = ((float) std::rand()) / RAND_MAX;
            int error_type;
            if (random <= mismatch_rate) {
                error_type = 0;
            } else {
                error_type = rand_iid(1,3);
            }
            switch (error_type) {
                case 0:
                    generate_candidate_text_add_mismatch(candidate_text,candidate_length);
                    break;
                case 1:
                    generate_candidate_text_add_deletion(candidate_text,&candidate_length);
                    break;
                default:
                    generate_candidate_text_add_insertion(candidate_text,&candidate_length);
                    break;
            }
        }
        // Return candidate-text
        candidate_text[candidate_length] = '\0';
        return candidate_text;
    }

public:
    Dataset(int _num_reads, int _length, float _error_rate, float _mismatch_rate, bool _exact_error_rate = false, bool _overwrite = false) {
        std::srand(static_cast<unsigned int>(std::time(nullptr)));
        num_reads = _num_reads;
        length = _length;
        error_rate = _error_rate;
        if (error_rate < 0 || error_rate > 0.7) {
            fprintf(stderr, "The error rate %f is either too small or too large. It should be between 0 and 0.7",
                    error_rate);
            exit(1);
        }
        mismatch_rate = _mismatch_rate;
        if (mismatch_rate < 0 || mismatch_rate > 1) {
            fprintf(stderr, "The mismatch rate %f is either too small or too large. It should be between 0 and 1",
                    mismatch_rate);
            exit(1);
        }
        exact_error_rate = _exact_error_rate;
        overwrite = _overwrite;
    }

    void output(const char * output_dir) {
        FILE *output_file;
        // Open files
        if (!overwrite) {
            if (FILE *file = fopen(output_dir, "r")) {
                fclose(file);
                return;
            }

        }
        output_file = fopen(output_dir,"w");
        // Allocate
        char* const pattern = new char[length + 1];
        // Read-align loop
        srand(time(0));
        int i;
        for (i = 0; i < num_reads; ++i) {
            // Prepare pattern
            generate_pattern(pattern, length);
            // Print pattern
            fprintf(output_file, ">%s\n", pattern);
            // Generate candidate-text
            char* const candidate_text = generate_candidate_text(
                    pattern, length, error_rate);
            // Print candidate-text
            fprintf(output_file, "<%s\n", candidate_text);
            delete[] candidate_text;
        }
        // Close files & free
        fclose(output_file);
        delete[] pattern;
    }

    std::string output() {
        std::string output_file_name;
        output_file_name += "simulated_" + std::to_string(num_reads) + "_" + std::to_string(length) + "_" +
                            std::to_string(error_rate) + "_";
        if (exact_error_rate) {
            output_file_name += "eq.seq";
        } else {
            output_file_name += "lt_eq.seq";
        }
        output(output_file_name.c_str());
        return output_file_name;
    }
};

#endif //GASMA_BENCHMARK_DATASET_H
