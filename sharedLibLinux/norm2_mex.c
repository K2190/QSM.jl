#include <inttypes.h>
#include <math.h>
#include "norm2_mex.h"
#include <stddef.h>
#define BLOCK_SIZE 128


float
norm2f(const float *x, const uint8_t *G, const size_t n)
{
    return sqrtf(norm22f(x, G, n));
}


float
norm22f(const float *x, const uint8_t *G, const size_t n)
{
    float s;

    if (n <= BLOCK_SIZE) {
        s = G[0] ? x[0]*x[0] : 0.0f;
        for(size_t i = 1; i < n; ++i) {
            if (G[i]) {
                s += x[i]*x[i];
            }
        }

    } else {
        size_t m = n >> 1;
        s = norm22f(x, G, m) + norm22f(x+m, G+m, n-m);
    }

    return s;
}


double
norm2d(const double *x, const uint8_t *G, const size_t n)
{
    return sqrt(norm22d(x, G, n));
}


double
norm22d(const double *x, const uint8_t *G, const size_t n)
{
    double s;

    if (n <= BLOCK_SIZE) {
        s = G[0] ? x[0]*x[0] : 0.0;
        for(size_t i = 1; i < n; ++i) {
            if (G[i]) {
                s += x[i]*x[i];
            }
        }

    } else {
        size_t m = n >> 1;
        s = norm22d(x, G, m) + norm22d(x+m, G+m, n-m);
    }

    return s;
}
