cmake_minimum_required(VERSION 3.16)
project(LEAP_SIMD)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mbmi -mavx2 -msse4.2")


set(LEAP_FILES
        print.cpp
        print.h
        shift.cpp
        shift.h
        popcount.cpp
        popcount.h
        SIMD_ED.cpp
        SIMD_ED.h
        SHD.cpp
        SHD.h
        bit_convert.cpp
        bit_convert.h
        mask.cpp
        mask.h)

add_executable(LEAP_SIMD main.cpp)
add_library(LEAP ${LEAP_FILES})
target_link_libraries(LEAP_SIMD LEAP)