#include <inttypes.h>
#include <omp.h>
#include <stddef.h>

void
gradfmf(float *dx, float *dy, float *dz,
        const float *u, const uint8_t *G,
        const double *h, const size_t *sz)
{
    size_t i, j, k;
    size_t l;

    const size_t nx = sz[0];
    const size_t ny = sz[1];
    const size_t nz = sz[2];
    const size_t nxny = nx*ny;
    const size_t nxnynz = nx*ny*nz;

    const size_t NX = nx-1;
    const size_t NY = nx*(ny-1);
    const size_t NZ = nxny*(nz-1);

    const float hx = (float)(1.0/h[0]);
    const float hy = (float)(1.0/h[1]);
    const float hz = (float)(1.0/h[2]);

    #pragma omp parallel for private(i,j,k,l) schedule(static) \
        if(nxnynz > 16*16*16)
    for(k = 0; k < nxnynz; k += nxny) {
        for(j = 0; j < nxny; j += nx) {
            l = j + k;
            for(i = 0; i < nx; ++i, ++l) {
                if (G[l]) {
                    dx[l] =
                        (i < NX) && G[l+1] ? hx*(u[l+1]-u[l]) :
                        (i > 0) && G[l-1] ? hx*(u[l]-u[l-1]) :
                        0.0f;

                    dy[l] =
                        (j < NY) && G[l+nx] ? hy*(u[l+nx]-u[l]) :
                        (j > 0) && G[l-nx] ? hy*(u[l]-u[l-nx]) :
                        0.0f;

                    dz[l] =
                        (k < NZ) && G[l+nxny] ? hz*(u[l+nxny]-u[l]) :
                        (k > 0) && G[l-nxny] ? hz*(u[l]-u[l-nxny]) :
                        0.0f;
                }
            }
        }
    }

    return;
}


void
gradfmd(double *dx, double *dy, double *dz,
        const double *u, const uint8_t *G,
        const double *h, const size_t *sz)
{
    size_t i, j, k;
    size_t l;

    const size_t nx = sz[0];
    const size_t ny = sz[1];
    const size_t nz = sz[2];
    const size_t nxny = nx*ny;
    const size_t nxnynz = nx*ny*nz;

    const size_t NX = nx-1;
    const size_t NY = nx*(ny-1);
    const size_t NZ = nxny*(nz-1);

    const double hx = 1.0/h[0];
    const double hy = 1.0/h[1];
    const double hz = 1.0/h[2];

    #pragma omp parallel for private(i,j,k,l) schedule(static) \
        if(nxnynz > 16*16*16)
    for(k = 0; k < nxnynz; k += nxny) {
        for(j = 0; j < nxny; j += nx) {
            l = j + k;
            for(i = 0; i < nx; ++i, ++l) {
                if (G[l]) {
                    dx[l] =
                        (i < NX) && G[l+1] ? hx*(u[l+1]-u[l]) :
                        (i > 0) && G[l-1] ? hx*(u[l]-u[l-1]) :
                        0.0;

                    dy[l] =
                        (j < NY) && G[l+nx] ? hy*(u[l+nx]-u[l]) :
                        (j > 0) && G[l-nx] ? hy*(u[l]-u[l-nx]) :
                        0.0;

                    dz[l] =
                        (k < NZ) && G[l+nxny] ? hz*(u[l+nxny]-u[l]) :
                        (k > 0) && G[l-nxny] ? hz*(u[l]-u[l-nxny]) :
                        0.0;
                }
            }
        }
    }

    return;
}
