#include <iostream>
#include "utils.h"

int main() {
    printf("wow");
    auto* matrix = new hurdle_matrix(
            "AGAGCTAAACATGGCCGCACATAAATCGTTTTGAGTTGAAACTTTACCGCTGCATCTATTTTTCTCCTAGAATTATACCGTACACAGCCGACGTTCCACC",
            "AGAGCTAAACAAGGGGCCCACATTAACGTTTTGAGCTTGAAGATCTTTACCGCGATCTATTTTTTCTCCTAGATTACCGTACACACCGACACTTCCATC",
            2, 128
            );
    matrix->print();
    return 0;
}
