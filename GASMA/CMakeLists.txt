# cmake_minimum_required(VERSION <specify CMake version here>)
project(GASMA)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse4.1 -DDEBUG")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1 -DDEBUG")

add_executable(GASMA main.cpp utils.h bit_convert.h bit_convert.cpp)