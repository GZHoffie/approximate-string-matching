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



add_executable(hurdle-matrix test/test_hurdle_matrix.cpp ${SHARED_FILES} hurdle_matrix.h)
add_executable(benchmark benchmark/benchmark.cpp  ${SHARED_FILES} hurdle_matrix.h)

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDISPLAY")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDISPLAY")

add_executable(test main.cpp utils.h bit_convert.h bit_convert.cpp mask.cpp mask.h hurdle_matrix.h)
