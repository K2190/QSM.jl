#pragma once

#include <stddef.h>

extern float norm2f(const float* x, const uint8_t* G, const size_t n);
extern double norm2d(const double* x, const uint8_t* G, const size_t n);
extern float norm22f(const float* x, const uint8_t* G, const size_t n);
extern double norm22d(const double* x, const uint8_t* G, const size_t n);