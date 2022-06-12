#include <inttypes.h>
#include <omp.h>
#include <stddef.h>
#include "lapmg_mex.h"
#include <stdio.h>


void
lapmgf(float *du,
       const float *u, const uint8_t *G, const double *h, const size_t *sz)
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

    #pragma omp parallel for private(i,j,k,l) schedule(static) \
        if(nxny*nz > 32*32*32)
    for(k = nxny; k < NZ; k += nxny) {
        for(j = nx; j < NY; j += nx) {
            l = 1 + j + k;
            for(i = 1; i < NX; ++i, ++l) {
                if (G[l]) {
                    du[l] =
                        hh*u[l] +
                        hx*(u[l-1] + u[l+1]) +
                        hy*(u[l-nx] + u[l+nx]) +
                        hz*(u[l-nxny] + u[l+nxny]);
                }
            }
        }
    }

    return;
}


void
lapmgd(double *du,
       const double *u, const uint8_t *G, const double *h, const size_t *sz)
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
    
    // printf("INSIDE lapmgd: nx= %d  ny= %d  nz= %d\n", nx, ny, nz);
    // printf("INSIDE lapmgd: hx= %f  hy= %f  hz= %f\n", hx, hy, hz);

    #pragma omp parallel for private(i,j,k,l) schedule(static) \
        if(nxny*nz > 32*32*32)
    for(k = nxny; k < NZ; k += nxny) {
        for(j = nx; j < NY; j += nx) {
            l = 1 + j + k;
            for(i = 1; i < NX; ++i, ++l) {
                if (G[l]) {
                    du[l] =
                        hh*u[l] +
                        hx*(u[l-1] + u[l+1]) +
                        hy*(u[l-nx] + u[l+nx]) +
                        hz*(u[l-nxny] + u[l+nxny]);
                }
            }
        }
    }

    return;
}
