# cmake_minimum_required(VERSION <specify CMake version here>)
project(GASMA)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mbmi -m64 -g3 -Wall -msse4.2 -march=native -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi -m64 -g3 -Wall -msse4.2 -march=native -O3")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")

# Options
option(DEBUG "display all the procedures of the algorithm" OFF)
option(DISPLAY "display the aligned strings" OFF)


add_executable(GASMA-hurdle-matrix test/test_hurdle_matrix.cpp utils.h bit_convert.h bit_convert.cpp mask.cpp mask.h hurdle_matrix.h hurdle_matrix_flipped.h)

add_executable(GASMA-hurdle-matrix-flipped test/test_hurdle_matrix_flipped.cpp utils.h bit_convert.h bit_convert.cpp mask.cpp mask.h hurdle_matrix.h hurdle_matrix_flipped.h)


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDISPLAY")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDISPLAY")

add_executable(GASMA main.cpp utils.h bit_convert.h bit_convert.cpp mask.cpp mask.h hurdle_matrix.h hurdle_matrix_flipped.h)
