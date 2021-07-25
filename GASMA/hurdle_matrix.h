//
// Created by 顾振昊 on 2021/5/22.
//

#ifndef APPROXIMATE_STRING_MATCHING_HURDLE_MATRIX_H
#define APPROXIMATE_STRING_MATCHING_HURDLE_MATRIX_H

#include <x86intrin.h>
#include <string>

/**
 * The maximum length of each string to be processed is 256
 * If the strings are longer than 256, then we need to use several
 * strings of length <= 256 to store it.
 */
#define _MAX_LENGTH 256

/**
 * The hurdle_matrix class
 * Implements a way to store the hurdle information of lanes using an
 * efficient way.
 */
class hurdle_matrix {
private:
    int k;                       // The k of banded Levenshtein distance
    bool remove_single_zeros;    // Whether we remove single zeros inside matrix
    bool initialized;            // Whether this matrix is initialized

    __m256i *mat;                // The matrix storing hurdle information

    /**
     * Load two DNAs, read and ref, into the class and calculate the
     * resulting hurdle matrix.
     *
     * @param read a string containing only 'A', 'C', 'G' and 'T'
     * @param ref a string containing only 'A', 'C', 'G' and 'T'
     * @param shift an integer indicating which lane in the hurdle matrix
     *     are we considering
    */
    __m256i hurdle_info(const std::string &read, const std::string &ref, int shift);


public:
    /**
     * Constructor of hurdle_matrix class
     * 
     * @param k the k of banded Levenshtein problem.
     * @param remove_single_zeros whether we remove the single zeros
     *     inside the hurdle matrix
     */
    hurdle_matrix(int k, bool remove_single_zeros);

    hurdle_matrix(const hurdle_matrix &) = default;
    
    ~hurdle_matrix();

    /** 
     * Load two DNAs, read and ref, into the class and calculate the
     * resulting hurdle matrix.
     * 
     * @param read a string containing only 'A', 'C', 'G' and 'T'
     * @param ref a string containing only 'A', 'C', 'G' and 'T'
    */
    void load_reads(const std::string &read, const std::string &ref);





};


#endif //APPROXIMATE_STRING_MATCHING_HURDLE_MATRIX_H
