#include <iostream>
#include "hurdle_matrix.h"

int main() {
    auto* matrix = new hurdle_matrix(
    "GCCCCGTCCCAGCACAGGCAGCGGGGATGCTCTAGAGCGATGTCGACCTGGGAAAGGCGCTGGGGCGCGCCGAATTCCAAAGGAGTTCCGTAAGGTTCAG",
    "GCTCCTCCCAGCACAGGCAGCGGGATGCTCTAGAGCATTCGCCCTGGGAAGGCGTGGGGCCCCCTAATTCCAAAAGGAGTTGCCCGTAAGGTTCAG"
            ,2, GLOBAL, 1, 3, 1);
    matrix->print();
    matrix->run();
    std::cout << "\n" << matrix->get_CIGAR() << "\n";
    std::cout << "cost: " << matrix->get_cost() << "\n";

    return 0;
}
