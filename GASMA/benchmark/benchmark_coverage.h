//
// Created by Zhenhao on 3/1/2022.
//

#ifndef GASMA_BENCHMARK_COVERAGE_H
#define GASMA_BENCHMARK_COVERAGE_H

#include <sstream>
#include <string>

/**
 * Find the long consecutive matching substring (LCM) of s1 and s2. A
 * substring of s1 and s2 is considered as an LCM if
 * 1) each of the characters match in s1 and s2,
 * 2) the length of the substring is larger or equal to threshold.
 *
 * @requires CIGAR string is obtained from aligning s1 and s2.
 * Supported CIGAR characters: 'I' - insert, 'D' - delete, 'X' - mismatch, '='/'M' - match.
 *
 * @param s1 the read string.
 * @param s2 the reference string.
 * @param CIGAR The CIGAR string.
 * @param threshold the smallest length of substring to be considered as LCM.
 * @return The LCM of s1 and s2.
 */
std::string long_consecutive_matching_substring(
        const char* s1,
        const char* s2,
        const std::string& CIGAR,
        int threshold = 3
        ) {
    std::string LCM;
    std::istringstream CIGAR_stream(CIGAR);
    int s1_index = 0;
    int s2_index = 0;

    int length;
    char type;
    while (CIGAR_stream >> length >> type) {
        switch (type) {
            case 'X':
                s1_index += length;
                s2_index += length;
                break;
            case 'I':
                s1_index += length;
                break;
            case 'D':
                s2_index += length;
                break;
            case '=':
            case 'M':
                for (int i = 0; i < length; i++) {
                    if (length >= threshold) {
                        LCM += s1[s1_index];
                    }
                    s1_index++;
                    s2_index++;

                }
                break;
            default:
                break;

        }

    }
    return LCM;


}


#endif //GASMA_BENCHMARK_COVERAGE_H
