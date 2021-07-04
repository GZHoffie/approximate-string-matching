# cmake_minimum_required(VERSION <specify CMake version here>)
project(GASMA)

set(CMAKE_CXX_STANDARD 14)

add_executable(GASMA main.cpp hurdle_matrix.h hurdle_matrix.cpp utils.h)