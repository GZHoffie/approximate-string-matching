//
// Created by Zhenhao on 11/5/2022.
//

#include <span>
#include "seqan3/alphabet/nucleotide/dna5.hpp"

/**
 * Convert a vector of seqan3::dna5 to string
 * @param dna vector of seqan3::dna5
 * @return the corresponding string
 */
std::string seqan_dna_to_cstring(std::span<seqan3::dna5> dna) {
    std::string res;
    res.reserve(dna.size());
    for (seqan3::dna5 c : dna) {
        res += c.to_char();
    }
    return res;
}