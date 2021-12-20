#include <iostream>
#include "hurdle_matrix.h"
#include "hurdle_matrix_flipped.h"

int main() {
    auto* matrix = new hurdle_matrix(
    "AGATAAACGAGGGTTGTAGTGGGACAATAGTCCAACACTTGCCACCTCCCAATGAATAAATACAAGCCTAAGGCGATCCATCCGACTTGATCAACCGGGG",
    "AGATATACGAGGTTGCGAGGGGGCAACTATTGTCCAAAACACTTAGCCACTATCCAATGATATAAACAAGCCTAGAGGCGATCATTGACTGATCACCGGG"
            ,3);
    matrix->print();
    matrix->run();
    std::cout << "\n" << matrix->get_CIGAR() << "\n";
    std::cout << "cost: " << matrix->get_cost() << "\n";

    return 0;
}
