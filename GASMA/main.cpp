#include <iostream>
#include "hurdle_matrix_flipped.h"

int main() {
    auto* matrix = new hurdle_matrix_flipped(
                                             "GCCCCGTCCCAGCACAGGCAGCGGGGATGCTCTAGAGCGATGTCGACCTGGGAAAGGCGCTGGGGCGCGCCGAATTCCAAAGGAGTTCCGTAAGGTTCAG",
                                             "GCTCCTCCCAGCACAGGCAGCGGGATGCTCTAGAGCATTCGCCCTGGGAAGGCGTGGGGCCCCCTAATTCCAAAAGGAGTTGCCCGTAAGGTTCAG",
                                             7);
    matrix->print();
    matrix->run();
    std::cout << "\n" << matrix->get_CIGAR() << "\n";
    delete matrix;
    return 0;
}
