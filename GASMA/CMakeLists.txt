cmake_minimum_required(VERSION 3.16)
project(GASMA)


set(CMAKE_CXX_STANDARD 14)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mbmi -m64 -g3 -Wall -msse4.2 -march=native -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi -m64 -g3 -Wall -msse4.2 -march=native -O3")

LINK_DIRECTORIES(
        benchmark/parasail/build
        benchmark/parasail/parasail
        benchmark/parasail
        benchmark/LEAP
)

# set file sets
set(SHARED_FILES
        ./utils.h
        ./bit_convert.h
        ./bit_convert.cpp
        ./mask.cpp ./mask.h)


# Add parasail library
ADD_SUBDIRECTORY(benchmark/parasail)

# Add LEAP library
ADD_SUBDIRECTORY(benchmark/LEAP_SIMD)


add_executable(hurdle-matrix test/test_hurdle_matrix.cpp ${SHARED_FILES} hurdle_matrix.h)
add_executable(hurdle-matrix-benchmark benchmark/benchmark.cpp ${SHARED_FILES} hurdle_matrix.h )

TARGET_LINK_DIRECTORIES(hurdle-matrix-benchmark PUBLIC
        benchmark/parasail/build
        benchmark/parasail/parasail
        benchmark/parasail
        benchmark/LEAP_SIMD
)

TARGET_LINK_LIBRARIES(hurdle-matrix-benchmark LEAP parasail)




#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDISPLAY")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDISPLAY")

add_executable(test main.cpp utils.h bit_convert.h bit_convert.cpp mask.cpp mask.h hurdle_matrix.h)
