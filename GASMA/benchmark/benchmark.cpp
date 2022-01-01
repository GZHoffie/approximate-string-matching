//
// Created by Zhenhao on 31/12/2021.
//

#include "benchmark_utils.h"

int main () {
    auto bench = new benchmark(1, 1, 1, 3, 100000, true);
    bench->read_string_file("../../pymatch/test/resource/sample.random.dataset.seq");
    bench->read_answer_file("../../pymatch/test/resource/nw2.txt");
    bench->run();
    bench->print();
}