# cmake_minimum_required(VERSION <specify CMake version here>)
project(GASMA)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse4.1 -march=native -O3 -DDISPLAY")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1 -march=native -O3 -DDISPLAY")

add_executable(GASMA-hurdle-matrix test/test_hurdle_matrix.cpp utils.h bit_convert.h bit_convert.cpp mask.cpp mask.h hurdle_matrix.h hurdle_matrix_flipped.h)

add_executable(GASMA-hurdle-matrix-flipped test/test_hurdle_matrix_flipped.cpp utils.h bit_convert.h bit_convert.cpp mask.cpp mask.h hurdle_matrix.h hurdle_matrix_flipped.h)


#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDEBUG")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEBUG")

add_executable(GASMA main.cpp utils.h bit_convert.h bit_convert.cpp mask.cpp mask.h hurdle_matrix.h hurdle_matrix_flipped.h)
