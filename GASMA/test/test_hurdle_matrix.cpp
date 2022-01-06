#include <iostream>
#include <string>
#include <sys/times.h>
#include <unistd.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include "../hurdle_matrix.h"

#define STRING_DIR "../../pymatch/test/resource/sample.random.dataset.seq"
#define ANSWER_DIR "../../pymatch/test/resource/nw2.txt"
#define TEST_NUM 100000
#define LOWER_ERROR_LIMIT 0
#define UPPER_ERROR_LIMIT 20

int main() {
    std::ifstream string_file, answer_file;
    auto *read = new std::string[TEST_NUM];
    auto *ref = new std::string[TEST_NUM];
    auto *optimal_res = new int[TEST_NUM];
    std::string read_temp, ref_temp;
    std::string optimal_res_temp;
    int pass = 0;
    int total = 0;
    string_file.open(STRING_DIR);
    answer_file.open(ANSWER_DIR);

    tms start_time;
    tms end_time;
    tms elp_time;

    elp_time.tms_stime = 0;
    elp_time.tms_utime = 0;
    elp_time.tms_cstime = 0;
    elp_time.tms_cutime = 0;

    auto* matrix = new hurdle_matrix;

    for (int i = 0; i < TEST_NUM; i++) {
        string_file.ignore();
        if (!getline(string_file, read_temp)) {
            break;
        }
        string_file.ignore();
        getline(string_file, ref_temp);
        read[i].assign(read_temp);
        ref[i].assign(ref_temp);

        getline(answer_file, optimal_res_temp);
        optimal_res[i] = atoi(optimal_res_temp.c_str());
        if (optimal_res[i] <= UPPER_ERROR_LIMIT && optimal_res[i] >= LOWER_ERROR_LIMIT) {
            total += 1;
        }
    }

    string_file.close();
    answer_file.close();


    for (int i = 0; i < TEST_NUM; i++) {
        //std::cout << read[i].c_str() << std::endl;
        //std::cout << ref[i].c_str() << std::endl;
        const char* s1 = read[i].c_str();
        const char* s2 = ref[i].c_str();
        int s1Len = (int) read[i].length();
        int s2Len = (int) ref[i].length();
        times(&start_time);
        matrix->reset(s1, s1Len, s2, s2Len, 2);
        matrix->run();
        times(&end_time);
        printf("%d %d\n", matrix->get_cost(), optimal_res[i]);
        if (optimal_res[i] <= UPPER_ERROR_LIMIT && optimal_res[i] >= LOWER_ERROR_LIMIT) {
            elp_time.tms_stime += end_time.tms_stime - start_time.tms_stime;
            elp_time.tms_utime += end_time.tms_utime - start_time.tms_utime;
            if (matrix->get_cost() == optimal_res[i]) {
                pass += 1;
            }
        }

    }



    printf("avg time: %f\n", (double) elp_time.tms_utime * TEST_NUM / total / sysconf(_SC_CLK_TCK));
    printf("total num: %d\n", total);
    printf("pass num: %d\n", pass);
    printf("correct rate: %f\n", (double) pass / (double) total);

    delete matrix;
    delete[] read;
    delete[] ref;
    delete[] optimal_res;
    return 0;
}
