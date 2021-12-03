#include <iostream>
#include "utils.h"

int main() {
    auto* matrix = new hurdle_matrix(
            "AGAGCTAAACATGGCCGCACATAAATCGTTTTGAGTTGAAACTTTACCGCTGCATCTATTTTTCTCCTAGAATTATACCGTACACAGCCGACGTTCCACC",
            "AGAGCTAAACAAGGGGCCCACATTAACGTTTTGAGCTTGAAGATCTTTACCGCGATCTATTTTTTCTCCTAGATTACCGTACACACCGACACTTCCATC",
            3, 128
            );
    matrix->print();
    matrix->run();
    delete matrix;
    return 0;
}
