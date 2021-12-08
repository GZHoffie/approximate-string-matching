#include <iostream>
#include "utils.h"

int main() {
    auto* matrix = new hurdle_matrix(
            "GCCCCGTCCCAGCACAGGCAGCGGGGATGCTCTAGAGCGATGTCGACCTGGGAAAGGCGCTGGGGCGCGCCGAATTCCAAAGGAGTTCCGTAAGGTTCAG",
            "GCTCCTCCCAGCACAGGCAGCGGGATGCTCTAGAGCATTCGCCCTGGGAAGGCGTGGGGCCCCCTAATTCCAAAAGGAGTTGCCCGTAAGGTTCAG",
            7
            );
    matrix->print();
    matrix->run();
    delete matrix;
    return 0;
}
