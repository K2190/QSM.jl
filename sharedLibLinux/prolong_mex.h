#pragma once

#include <stddef.h>

extern void prolongf(float* x,
    const float* x2, const uint8_t* G,
    const size_t* sz, const size_t* sz2);

extern void prolongd(double* x,
    const double* x2, const uint8_t* G,
    const size_t* sz, const size_t* sz2);

