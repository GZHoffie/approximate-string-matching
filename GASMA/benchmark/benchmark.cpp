//
// Created by Zhenhao on 31/12/2021.
//

#include "benchmark_utils.h"
#include "benchmark_dataset.h"

#include <vector>

int main () {
    int num_reads = 5000000;
    int length = 100;
    std::vector<float> error_rates = {0.05, 0.10, 0.15, 0.20};

    for (auto error_rate : error_rates) {
        Dataset dataset(num_reads, length, error_rate, true);
        std::string output_dir = dataset.output();

        benchmark bench(1, 1, 1, 10, 1000000, true );
        bench.read_string_file(output_dir.c_str());
        bench.run();
        bench.print();
    }

}