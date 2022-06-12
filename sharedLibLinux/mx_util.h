#pragma once

#include <inttypes.h>
#include <stdint.h>

extern void     mx_zero(double* mxx);
double* mx_pad_boundary(const size_t* sz, const double* mxx, size_t* new_sz);

uint8_t* mx_pad_boundary_uint(const size_t* sz, const uint8_t* mxx, size_t* new_sz);


extern double* mx_unpad_boundary(const double* mxxp, size_t* mxxp_sz, double* mxx);
