#include <iostream>
#include "hurdle_matrix.h"
#include "benchmark/benchmark_coverage.h"

int main() {
    auto* matrix = new hurdle_matrix<int_128bit>(
    "TGCGTCGCAGTGTCGAAAAAAATCTAAAGGGCCTACCGGGCCCCTTCATGAAGAGCAGGGAAAGAATCTCTATACTATTACTAATGGAAAGCAGTCTTGG",
    "TGCGTCGCAGTGTCGAAATAAAATCTAAAGGGCCTACCGGGCCCCTTCATGAAGAGTAGGTAAAGAATCTCTAGTACTACTACTAATGGAAAGCAGTCTTGG"
            ,3, GLOBAL, 1, 1, 1);
    matrix->print();
    matrix->run();
    std::cout << "\n" << matrix->get_CIGAR() << "\n";
    std::cout << "cost: " << matrix->get_cost() << "\n";
    std::cout << long_consecutive_matching_substring("GCCCCGTCCCAGCACAGGCAGCGGGGATGCTCTAGAGCGATGTCGACCTGGGAAAGGCGCTGGGGCGCGCCGAATTCCAAAGGAGTTCCGTAAGGTTCAG",
                                                     "GCTCCTCCCAGCACAGGCAGCGGGATGCTCTAGAGCATTCGCCCTGGGAAGGCGTGGGGCCCCCTAATTCCAAAAGGAGTTGCCCGTAAGGTTCAG",
                                                     matrix->get_CIGAR()) << std::endl;

    return 0;
}
