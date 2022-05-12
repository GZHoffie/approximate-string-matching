cmake_minimum_required(VERSION 3.4)
PROJECT(GASMA CXX)

SET(CMAKE_CXX_STANDARD 20)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mbmi -m64 -g3 -Wall -msse4.2 -march=native -O3")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi -m64 -g3 -Wall -msse4.2 -march=native -O3")

LINK_DIRECTORIES(
        benchmark/parasail/build
        benchmark/parasail/parasail
        benchmark/parasail
        benchmark/LEAP
        mapper/seqan3/submodules/cereal
)

FIND_PACKAGE(SeqAn3 3.0.0 REQUIRED HINTS "${CMAKE_SOURCE_DIR}/mapper/seqan3/build_system")

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

# Add seqan3 library
ADD_SUBDIRECTORY(mapper/seqan3)
ADD_SUBDIRECTORY(mapper/seqan3/submodules/cereal)

# Executable for running hurdle-matrix
ADD_EXECUTABLE(hurdle-matrix test/test_hurdle_matrix.cpp ${SHARED_FILES} hurdle_matrix.h)
SET_TARGET_PROPERTIES(hurdle-matrix PROPERTIES COMPILE_FLAGS "-DDISPLAY")

# Executable for Benchmarking
ADD_EXECUTABLE(hurdle-matrix-benchmark benchmark/benchmark.cpp ${SHARED_FILES} hurdle_matrix.h benchmark/benchmark_coverage.h benchmark/benchmark_dataset.h)
#SET_TARGET_PROPERTIES(hurdle-matrix-benchmark PROPERTIES COMPILE_FLAGS "-DDEBUG -DDISPLAY")
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
ADD_LIBRARY(GASMA ${SHARED_FILES} hurdle_matrix.h main.cpp)

# Executable for mapper
ADD_EXECUTABLE(my-mapper ${SHARED_FILES} mapper/main.cpp seqan3_main.h)
TARGET_LINK_LIBRARIES(my-mapper seqan3::seqan3 cereal)

ADD_EXECUTABLE(my-indexer mapper/indexer.cpp)
TARGET_LINK_LIBRARIES(my-indexer seqan3::seqan3)