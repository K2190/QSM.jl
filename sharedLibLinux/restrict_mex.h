#pragma once

#include <stddef.h>
#include <inttypes.h>

extern void restrictf(float* x2,
    const float* x, const uint8_t* G2,
    const size_t* sz2, const size_t* sz);

extern void restrictd(double* x2,
    const double* x, const uint8_t* G2,
    const size_t* sz2, const size_t* sz);
