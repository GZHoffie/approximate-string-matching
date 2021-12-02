#include <iostream>
#include "utils.h"

int main() {
    auto* matrix = new hurdle_matrix(
            "AGAGCTAAACATGGCCGCACATAAATCGTTTTGAGTTGAAACTTTACCGCTGCATCTATTTTTCTCCTAGAATTATACCGTACACAGCCGACGTTCCACC",
            "AGAGCTAAACAAGGGGCCCACATTAACGTTTTGAGCTTGAAGATCTTTACCGCGATCTATTTTTTCTCCTAGATTACCGTACACACCGACACTTCCATC",
            2, 128
            );
    matrix->print();
    delete matrix;
    return 0;
}
