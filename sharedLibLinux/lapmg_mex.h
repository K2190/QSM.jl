#pragma once



extern void
lapmgf(float *du,
       const float *u, const uint8_t *G, const double *h, const size_t *sz);

extern void
lapmgd(double *du,
       const double *u, const uint8_t *G, const double *h, const size_t *sz);
