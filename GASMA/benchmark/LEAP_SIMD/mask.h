//
// Created by Zhenhao on 2/1/2022.
//

#ifndef LEAP_SIMD_MASK_H
#define LEAP_SIMD_MASK_H

#include <cstdint>
/*
#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/iteration.hpp>
#include <boost/preprocessor/arithmetic.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>
*/

#ifndef __aligned
#define __aligned __attribute__((aligned(32)))
#endif

#define SSE_BIT_LENGTH		128
#define SSE_BYTE_LENGTH		16
#define AVX_BIT_LENGTH		256
#define AVX_BYTE_LENGTH		32
#define MAX_ERROR_SSE		8
#define MAX_ERROR_AVX		16

extern uint8_t *MASK_SSE_END;
extern uint8_t *MASK_SSE_BEG;
extern uint8_t *MASK_AVX_END;
extern uint8_t *MASK_AVX_BEG;

extern uint8_t *MASK_SRS;

#endif //LEAP_SIMD_MASK_H
