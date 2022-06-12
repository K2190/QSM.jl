#include <inttypes.h>
#include <omp.h>
#include <stddef.h>
#include "residual_mex.h"

void
residualf(float *r,
          const float *f, const float *x, const uint8_t *G,
          const double *h, const size_t *sz)
{
    size_t i, j, k;
    size_t l;

    const size_t nx = sz[0];
    const size_t ny = sz[1];
    const size_t nz = sz[2];
    const size_t nxny = nx*ny;

    const size_t NX = nx-1;
    const size_t NY = nx*(ny-1);
    const size_t NZ = nxny*(nz-1);

    const float hx = (float)(1.0/(h[0]*h[0]));
    const float hy = (float)(1.0/(h[1]*h[1]));
    const float hz = (float)(1.0/(h[2]*h[2]));
    const float hh = (float)(-2.0*(hx+hy+hz));

    // printf("INSIDE residualf :nx, ny, nz = : [%d] [%d] [%d]\n", nx, ny, nz);
    // printf("INSIDE residualf :hx, hy, hz = : [%f] [%f] [%f]\n", hx, hy, hz);

    #pragma omp parallel for private(i,j,k,l) schedule(static) \
        if (nxny*nz > 32*32*32)
    for(k = nxny; k < NZ; k += nxny) {
        for(j = nx; j < NY; j += nx) {
            l = 1 + j + k;
            for(i = 1; i < NX; ++i, ++l) {
                if (G[l]) {
                    r[l] = f[l] +
                        (hh*x[l] +
                        hx*(x[l-1] + x[l+1]) +
                        hy*(x[l-nx] + x[l+nx]) +
                        hz*(x[l-nxny] + x[l+nxny]));
                }
            }
        }
    }

    return;
}


void
residuald(double *r,
          const double *f, const double *x, const uint8_t *G,
          const double *h, const size_t *sz)
{
    size_t i, j, k;
    size_t l;

    const size_t nx = sz[0];
    const size_t ny = sz[1];
    const size_t nz = sz[2];
    const size_t nxny = nx*ny;

    const size_t NX = nx-1;
    const size_t NY = nx*(ny-1);
    const size_t NZ = nxny*(nz-1);

    const double hx = 1.0/(h[0]*h[0]);
    const double hy = 1.0/(h[1]*h[1]);
    const double hz = 1.0/(h[2]*h[2]);
    const double hh = -2.0*(hx+hy+hz);

    // printf("INSIDE residuald :nx, ny, nz = : [%d] [%d] [%d]\n", nx, ny, nz);
    // printf("INSIDE residuald :hx, hy, hz = : [%f] [%f] [%f]\n", hx, hy, hz);

    #pragma omp parallel for private(i,j,k,l) schedule(static) \
        if (nxny*nz > 32*32*32)
    for(k = nxny; k < NZ; k += nxny) {
        for(j = nx; j < NY; j += nx) {
            l = 1 + j + k;
            for(i = 1; i < NX; ++i, ++l) {
                if (G[l]) {
                    r[l] = f[l] +
                        (hh*x[l] +
                        hx*(x[l-1] + x[l+1]) +
                        hy*(x[l-nx] + x[l+nx]) +
                        hz*(x[l-nxny] + x[l+nxny]));
                }
            }
        }
    }

    return;
}
