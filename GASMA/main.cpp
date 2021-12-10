#include <iostream>
#include "utils.h"

int main() {
    auto* matrix = new hurdle_matrix(
    "CCGGTTACGATGGCTGCCGGAATCGGGGCTTTCCCACGGTTGCCACTCTGATGAACGTTTTTTAACCGTTCGGTTGCTAACATTAAATAAGTTGAATCGG",
    "CCGGTTACGATGGCGGCCGGAATCGGGCGCTTTCTCACGGTTGCCACTCTGATGTAGGTTTTTTAACCGTTCGGTTGCTACATTAAAATAAGTTGAATCGGG"
            ,2
            );
    matrix->print();
    matrix->run();
    std::cout << "\n" << matrix->get_CIGAR() << "\n";
    delete matrix;
    return 0;
}
