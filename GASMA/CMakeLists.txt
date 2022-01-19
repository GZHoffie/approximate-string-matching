cmake_minimum_required(VERSION 3.16)
PROJECT(GASMA)

SET(CMAKE_CXX_STANDARD 14)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mbmi -m64 -g3 -Wall -msse4.2 -march=native -O3")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi -m64 -g3 -Wall -msse4.2 -march=native -O3")

LINK_DIRECTORIES(
        benchmark/parasail/build
        benchmark/parasail/parasail
        benchmark/parasail
        benchmark/LEAP
)

# set file sets
SET(SHARED_FILES
        ./utils.h
        ./bit_convert.h
        ./bit_convert.cpp
        ./mask.cpp ./mask.h)


# Add parasail library
ADD_SUBDIRECTORY(benchmark/parasail)

# Add LEAP library
ADD_SUBDIRECTORY(benchmark/LEAP_SIMD)

# Executable for running hurdle-matrix
ADD_EXECUTABLE(hurdle-matrix test/test_hurdle_matrix.cpp ${SHARED_FILES} hurdle_matrix.h)
SET_TARGET_PROPERTIES(hurdle-matrix PROPERTIES COMPILE_FLAGS "-DDISPLAY")

# Executable for Benchmarking
ADD_EXECUTABLE(hurdle-matrix-benchmark benchmark/benchmark.cpp ${SHARED_FILES} hurdle_matrix.h benchmark/benchmark_coverage.h benchmark/benchmark_dataset.h)
TARGET_LINK_DIRECTORIES(hurdle-matrix-benchmark PUBLIC
        benchmark/parasail/build
        benchmark/parasail/parasail
        benchmark/parasail
        benchmark/LEAP_SIMD
)
TARGET_LINK_LIBRARIES(hurdle-matrix-benchmark LEAP parasail)

# Executable for testing functions
ADD_EXECUTABLE(test main.cpp utils.h bit_convert.h bit_convert.cpp mask.cpp mask.h hurdle_matrix.h benchmark/benchmark_coverage.h)
SET_TARGET_PROPERTIES(test PROPERTIES COMPILE_FLAGS "-DDEBUG -DDISPLAY")

# Compiling the library for greedy algorithm
ADD_LIBRARY(GASMA ${SHARED_FILES} hurdle_matrix.h)