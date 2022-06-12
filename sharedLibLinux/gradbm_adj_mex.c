#include <inttypes.h>
#include <omp.h>
#include <stddef.h>

void
gradbm_adjf(float *du,
            const float *x, const float *y, const float *z,
            const uint8_t *G, const double *h, const size_t *sz)
{
    size_t i, j, k;
    size_t l;

    float dx, dy, dz;

    const size_t nx = sz[0];
    const size_t ny = sz[1];
    const size_t nz = sz[2];
    const size_t nxny = nx*ny;
    const size_t nxnynz = nx*ny*nz;

    const size_t NX = nx-1;
    const size_t NY = nx*(ny-1);
    const size_t NZ = nxny*(nz-1);

    const float hx = (float)(-1.0/h[0]);
    const float hy = (float)(-1.0/h[1]);
    const float hz = (float)(-1.0/h[2]);

    #pragma omp parallel for private(i,j,k,l) schedule(static) \
        if(nxnynz > 16*16*16)
    for(k = 0; k < nxnynz; k += nxny) {
        for(j = 0; j < nxny; j += nx) {
            l = j + k;
            for(i = 0; i < nx; ++i, ++l) {
                if (G[l]) {
                    dx =
                        (i < NX) && G[l+1] ? hx*(x[l+1]-x[l]) :
                        (i > 0) && G[l-1] ? hx*(x[l]-x[l-1]) :
                        0.0f;

                    dy =
                        (j < NY) && G[l+nx] ? hy*(y[l+nx]-y[l]) :
                        (j > 0) && G[l-nx] ? hy*(y[l]-y[l-nx]) :
                        0.0f;

                    dz =
                        (k < NZ) && G[l+nxny] ? hz*(z[l+nxny]-z[l]) :
                        (k > 0) && G[l-nxny] ? hz*(z[l]-z[l-nxny]) :
                        0.0f;

                    du[l] = dx + dy + dz;
                }
            }
        }
    }

    return;
}


void
gradbm_adjd(double *du,
            const double *x, const double *y, const double *z,
            const uint8_t *G, const double *h, const size_t *sz)
{
    size_t i, j, k;
    size_t l;

    double dx, dy, dz;

    const size_t nx = sz[0];
    const size_t ny = sz[1];
    const size_t nz = sz[2];
    const size_t nxny = nx*ny;
    const size_t nxnynz = nx*ny*nz;

    const size_t NX = nx-1;
    const size_t NY = nx*(ny-1);
    const size_t NZ = nxny*(nz-1);

    const double hx = -1.0/h[0];
    const double hy = -1.0/h[1];
    const double hz = -1.0/h[2];

    #pragma omp parallel for private(i,j,k,l) schedule(static) \
        if(nxnynz > 16*16*16)
    for(k = 0; k < nxnynz; k += nxny) {
        for(j = 0; j < nxny; j += nx) {
            l = j + k;
            for(i = 0; i < nx; ++i, ++l) {
                if (G[l]) {
                    dx =
                        (i < NX) && G[l+1] ? hx*(x[l+1]-x[l]) :
                        (i > 0) && G[l-1] ? hx*(x[l]-x[l-1]) :
                        0.0;

                    dy =
                        (j < NY) && G[l+nx] ? hy*(y[l+nx]-y[l]) :
                        (j > 0) && G[l-nx] ? hy*(y[l]-y[l-nx]) :
                        0.0;

                    dz =
                        (k < NZ) && G[l+nxny] ? hz*(z[l+nxny]-z[l]) :
                        (k > 0) && G[l-nxny] ? hz*(z[l]-z[l-nxny]) :
                        0.0;

                    du[l] = dx + dy + dz;
                }
            }
        }
    }

    return;
}
