#pragma once

#include <inttypes.h>

extern void
gauss_seidelf(float* v,
    const float* f, const uint8_t* G,
    const double* h, const size_t* sz,
    int32_t iter, const uint8_t rev);

extern void
gauss_seideld(double* v,
    const double* f, const uint8_t* G,
    const double* h, int32_t iter, const uint8_t rev,
    const size_t* sz);

