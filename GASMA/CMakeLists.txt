# cmake_minimum_required(VERSION <specify CMake version here>)
project(GASMA)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mbmi -m64 -g3 -Wall -msse4.2 -march=native -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi -m64 -g3 -Wall -msse4.2 -march=native -O3")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")

# set file sets
set(SHARED_FILES
        utils.h
        bit_convert.h
        bit_convert.cpp
        mask.cpp mask.h)

set(BENCHMARK_SHARED_FILES benchmark/benchmark_utils.h)


add_executable(hurdle-matrix test/test_hurdle_matrix.cpp hurdle_matrix.h ${SHARED_FILES})
add_executable(hurdle-matrix-benchmark hurdle_matrix.h ${BENCHMARK_SHARED_FILES} ${SHARED_FILES})


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDISPLAY")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDISPLAY")

add_executable(test main.cpp utils.h bit_convert.h bit_convert.cpp mask.cpp mask.h hurdle_matrix.h)
