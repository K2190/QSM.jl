#pragma once

#include <stddef.h>

extern void residualf(float* r,
    const float* f, const float* x, const uint8_t* G,
    const double* h, const size_t* sz);

extern void residuald(double* r,
    const double* f, const double* x, const uint8_t* G,
    const double* h, const size_t* sz);

