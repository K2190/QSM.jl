#include <inttypes.h>
#include <omp.h>
#include <stddef.h>
#include "gauss_seidel_mex.h"
#include <stdio.h>

void
gauss_seidelf(float *v,
              const float *f, const uint8_t *G,
              const double *h, const size_t *sz,
              int32_t iter, const uint8_t rev)
{
    int32_t i, j, k;
    int32_t l;

    const int32_t nx = (int32_t)sz[0];
    const int32_t ny = (int32_t)sz[1];
    const int32_t nz = (int32_t)sz[2];
    const int32_t nxny = nx*ny;

    const int32_t NX = rev ? nx-2 : nx-1;
    const int32_t NY = rev ? nx*(ny-2) : nx*(ny-1);
    const int32_t NZ = rev ? nxny*(nz-2) : nxny*(nz-1);

    const float hx = (float)(1.0/(h[0]*h[0]));
    const float hy = (float)(1.0/(h[1]*h[1]));
    const float hz = (float)(1.0/(h[2]*h[2]));
    const float hh = (float)(1.0/(2.0*(double)(hx+hy+hz)));

#ifndef GAUSS_SEIDEL_RED_BLACK
    if (rev > 0) {
        while(iter-- > 0) {
            for(k = NZ; k >= nxny; k -= nxny) {
                for(j = NY; j >= nx; j -= nx) {
                    l = NX + j + k;
                    for(i = NX; i >= 1; --i, --l) {
                        if (G[l]) {
                            v[l] = hh *
                                (f[l] +
                                (hx*(v[l-1] + v[l+1]) +
                                hy*(v[l-nx] + v[l+nx]) +
                                hz*(v[l-nxny] + v[l+nxny])));
                        }
                    }
                }
            }
        }
    } else {
        while(iter-- > 0) {
            for(k = nxny; k < NZ; k += nxny) {
                for(j = nx; j < NY; j += nx) {
                    l = 1 + j + k;
                    for(i = 1; i < NX; ++i, ++l) {
                        if (G[l]) {
                            v[l] = hh *
                                (f[l] +
                                (hx*(v[l-1] + v[l+1]) +
                                hy*(v[l-nx] + v[l+nx]) +
                                hz*(v[l-nxny] + v[l+nxny])));
                        }
                    }
                }
            }
        }
    }

#else
    size_t s, is, js, ks;

    if (rev > 0) {
        while(iter-- > 0) {
            for(s = 0; s < 2; ++s) {
                #pragma omp parallel for private(i,j,k,l,is,js,ks) \
                    schedule(static) if (nxny*nz > 32*32*32)
                for(k = NZ; k >= nxny; k -= nxny) {
                    ks = (k/nxny) & 1;
                    for(j = NY; j >= nx; j -= nx) {
                        js = (j/nx) & 1;
                        is = s + ((ks && js) || !(js || ks));
                        l = NX-is + j + k;
                        for(i = NX-is; i >= 1; i -= 2, l -= 2) {
                            if (G[l]) {
                                v[l] = hh *
                                    (f[l] +
                                    (hx*(v[l-1] + v[l+1]) +
                                    hy*(v[l-nx] + v[l+nx]) +
                                    hz*(v[l-nxny] + v[l+nxny])));
                            }
                        }
                    }
                }
            }
        }
    } else {
        while(iter-- > 0) {
            for(s = 0; s < 2; ++s) {
                #pragma omp parallel for private(i,j,k,l,is,js,ks) \
                    schedule(static) if (nxny*nz > 32*32*32)
                for(k = nxny; k < NZ; k += nxny) {
                    ks = (k/nxny) & 1;
                    for(j = nx; j < NY; j += nx) {
                        js = (j/nx) & 1;
                        is = s + ((ks && js) || !(js || ks));
                        l = 1 + is + j + k;
                        for(i = 1 + is; i < NX; i += 2, l += 2) {
                            if (G[l]) {
                                v[l] = hh *
                                    (f[l] +
                                    (hx*(v[l-1] + v[l+1]) +
                                    hy*(v[l-nx] + v[l+nx]) +
                                    hz*(v[l-nxny] + v[l+nxny])));
                            }
                        }
                    }
                }
            }
        }
    }
#endif

    return;
}


void
gauss_seideld(double *v,
              const double *f, const uint8_t *G,
              const double *h, int32_t iter, const uint8_t rev,
              const size_t* sz)
{
    int32_t i, j, k;
    int32_t l;

    const int32_t nx = (int32_t)sz[0];
    const int32_t ny = (int32_t)sz[1];
    const int32_t nz = (int32_t)sz[2];
    const int32_t nxny = nx*ny;

    // printf("INSIDE gauss_seideld :nx, ny, nz = : [%d] [%d] [%d]\n", nx, ny, nz);

    const int32_t NX = rev ? nx-2 : nx-1;
    const int32_t NY = rev ? nx*(ny-2) : nx*(ny-1);
    const int32_t NZ = rev ? nxny*(nz-2) : nxny*(nz-1);

    const double hx = 1.0/(h[0]*h[0]);
    const double hy = 1.0/(h[1]*h[1]);
    const double hz = 1.0/(h[2]*h[2]);
    const double hh = 1.0/(2.0*(hx+hy+hz));

    // printf("INSIDE gauss_seideld :hx, hy, hz = : [%f] [%f] [%f]\n", hx, hy, hz);


#ifndef GAUSS_SEIDEL_RED_BLACK
    if (rev > 0) {
        while(iter-- > 0) {
            for(k = NZ; k >= nxny; k -= nxny) {
                for(j = NY; j >= nx; j -= nx) {
                    l = NX + j + k;
                    for(i = NX; i >= 1; --i, --l) {
                        if (G[l]) {
                            v[l] = hh *
                                (f[l] +
                                (hx*(v[l-1] + v[l+1]) +
                                hy*(v[l-nx] + v[l+nx]) +
                                hz*(v[l-nxny] + v[l+nxny])));
                        }
                    }
                }
            }
        }
    } else {
        while(iter-- > 0) {
            for(k = nxny; k < NZ; k += nxny) {
                for(j = nx; j < NY; j += nx) {
                    l = 1 + j + k;
                    for(i = 1; i < NX; ++i, ++l) {
                        if (G[l]) {
                            v[l] = hh *
                                (f[l] +
                                (hx*(v[l-1] + v[l+1]) +
                                hy*(v[l-nx] + v[l+nx]) +
                                hz*(v[l-nxny] + v[l+nxny])));
                        }
                    }
                }
            }
        }
    }

#else
    size_t s, is, js, ks;

    int sample_cnt = 0;
    
    if (rev > 0) {
        while(iter-- > 0) {
            for(s = 0; s < 2; ++s) {
                #pragma omp parallel for private(i,j,k,l,is,js,ks) \
                    schedule(static) if (nxny*nz > 32*32*32)
                for(k = NZ; k >= nxny; k -= nxny) {
                    ks = (k/nxny) & 1;
                    for(j = NY; j >= nx; j -= nx) {
                        js = (j/nx) & 1;
                        is = s + ((ks && js) || !(js || ks));
                        l = NX-is + j + k;
                        for(i = NX-is; i >= 1; i -= 2, l -= 2) {
                            if (G[l]) {
                                v[l] = hh *
                                    (f[l] +
                                    (hx*(v[l-1] + v[l+1]) +
                                    hy*(v[l-nx] + v[l+nx]) +
                                    hz*(v[l-nxny] + v[l+nxny])));

                                    if (0) { //l == 30276 || l == 7709) {
                                        printf("GaussSeidel DUMP 1\n");
                                        double sum = (hx*(v[l-1] + v[l+1]) + hy*(v[l-nx] + v[l+nx]) +  hz*(v[l-nxny] + v[l+nxny]));
                                        printf("SUM = %f\n", sum);
                                        printf(">>GaussSeidel inputs  hh = %f  hx = %f  hy = %f  hz = %f\n", hh, hx, hy, hz);
                                        printf(">>GaussSeidel inputs v[l-1] = %f   v[l+1] = %f\n", v[l-1], v[l+1]);
                                        printf(">>GaussSeidel inputs v[l-nx] = %f   v[l+nx] = %f\n", v[l-nx], v[l+nx]);
                                        printf(">>GaussSeidel inputs v[l-nxny] = %f   v[l+nxny] = %f\n", v[l-nxny], v[l+nxny]);
                                        printf(">>GaussSeidel inputs f[l] = %f   \n", f[l]);
                                        printf(">>GaussSeidel output p: l = %d G[l] = %d     v[l] = %f\n", l, G[l], v[l]);
                                        // sample_cnt++;
                                    }
                                    
                            }
                        }
                    }
                }
            }
        }
    } else {
        while(iter-- > 0) {
            for(s = 0; s < 2; ++s) {
                #pragma omp parallel for private(i,j,k,l,is,js,ks) \
                    schedule(static) if (nxny*nz > 32*32*32)
                for(k = nxny; k < NZ; k += nxny) {
                    ks = (k/nxny) & 1;
                    for(j = nx; j < NY; j += nx) {
                        js = (j/nx) & 1;
                        is = s + ((ks && js) || !(js || ks));
                        l = 1 + is + j + k;
                        for(i = 1 + is; i < NX; i += 2, l += 2) {
                            if (G[l]) {
                                v[l] = hh *
                                    (f[l] +
                                    (hx*(v[l-1] + v[l+1]) +
                                    hy*(v[l-nx] + v[l+nx]) +
                                    hz*(v[l-nxny] + v[l+nxny])));
                                if (0) { //l == 30276 || l == 7709) {
                                        double sum = (hx*(v[l-1] + v[l+1]) + hy*(v[l-nx] + v[l+nx]) +  hz*(v[l-nxny] + v[l+nxny]));
                                        printf("GaussSeidel DUMP 2\n");
                                        printf("SUM = %f\n", sum);
                                        printf(">>GaussSeidel inputs  hh = %f  hx = %f  hy = %f  hz = %f\n", hh, hx, hy, hz);
                                        printf(">>GaussSeidel inputs v[l-1] = %f   v[l+1] = %f\n", v[l-1], v[l+1]);
                                        printf(">>GaussSeidel inputs v[l-nx] = %f   v[l+nx] = %f\n", v[l-nx], v[l+nx]);
                                        printf(">>GaussSeidel inputs v[l-nxny] = %f   v[l+nxny] = %f\n", v[l-nxny], v[l+nxny]);
                                        printf(">>GaussSeidel inputs f[l] = %f   \n", f[l]);
                                        printf(">>GaussSeidel output p: l = %d G[l] = %d     v[l] = %f\n", l, G[l], v[l]);
                                       // sample_cnt++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
#endif

    return;
}
