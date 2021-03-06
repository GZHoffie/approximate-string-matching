/*
 * bit_convert.h
 *
 *  Created on: Nov 7, 2013
 *      Author: hxin
 */

#ifndef GASMA_BIT_CONVERT_H_
#define GASMA_BIT_CONVERT_H_

#include <cstdint>

#ifndef __aligned__
#define __aligned__ __attribute__((aligned(16)))
#endif

void c_convert2bit(char *str, int length, uint8_t *bits);

void sse3_convert2bit11(char *str, int length, uint8_t *bits);

void sse3_convert2bit1(char *str, uint8_t *bits0, uint8_t *bits1);

#endif /* BIT_CONVERT_H_ */