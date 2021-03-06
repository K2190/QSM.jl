#include <inttypes.h>
#include <omp.h>
#include <stddef.h>

void
gradf_adjf(float *du,
           const float *x, const float *y, const float *z,
           const double *h, const size_t *sz)
{
    size_t i, j, k;
    size_t l;

    const size_t nx = sz[0];
    const size_t ny = sz[1];
    const size_t nz = sz[2];
    const size_t nxny = nx*ny;
    const size_t nxnynz = nx*ny*nz;

    const float hx = (float)(-1.0/h[0]);
    const float hy = (float)(-1.0/h[1]);
    const float hz = (float)(-1.0/h[2]);

    /* i = 0, j = 0, k = 0 */
    l = 0;
    du[l] =
        hx*(x[l+1]-x[l]) +
        hy*(y[l+nx]-y[l]) +
        hz*(z[l+nxny]-z[l]);

#pragma omp parallel private(i,j,k,l) if (nxny*nz > 16*16*16)
{
    /* i = 0, k = 0 */
    #pragma omp for schedule(static)
    for(l = nx; l < nxny; l += nx) {
        du[l] =
            hx*(x[l+1]-x[l]) +
            hy*(y[l]-y[l-nx]) +
            hz*(z[l+nxny]-z[l]);
    }

    /* j = 0, k = 0 */
    #pragma omp for schedule(static)
    for(l = 1; l < nx; ++l) {
        du[l] =
            hx*(x[l]-x[l-1]) +
            hy*(y[l+nx]-y[l]) +
            hz*(z[l+nxny]-z[l]);
    }

    /* k = 0 */
    #pragma omp for schedule(static) collapse(2)
    for(j = nx; j < nxny; j += nx) {
        for(i = 1; i < nx; ++i) {
            l = i + j;
            du[l] =
                hx*(x[l]-x[l-1]) +
                hy*(y[l]-y[l-nx]) +
                hz*(z[l+nxny]-z[l]);
        }
    }

    /* interior loop */
    #pragma omp for schedule(static)
    for(k = nxny; k < nxnynz; k += nxny) {
        /* i = 0, j = 0 */
        l = k;
        du[l] =
            hx*(x[l+1]-x[l]) +
            hy*(y[l+nx]-y[l]) +
            hz*(z[l]-z[l-nxny]);

        /* j = 0 */
        l = 1 + k;
        for(i = 1; i < nx; ++i, ++l) {
            du[l] =
                hx*(x[l]-x[l-1]) +
                hy*(y[l+nx]-y[l]) +
                hz*(z[l]-z[l-nxny]);
        }

        for(j = nx; j < nxny; j += nx) {
            /* i = 0 */
            l = j + k;
            du[l] =
                hx*(x[l+1]-x[l]) +
                hy*(y[l]-y[l-nx]) +
                hz*(z[l]-z[l-nxny]);

            l = 1 + j + k;
            for(i = 1; i < nx; ++i, ++l) {
                du[l] =
                    hx*(x[l]-x[l-1]) +
                    hy*(y[l]-y[l-nx]) +
                    hz*(z[l]-z[l-nxny]);
            }
        }
    }

} /* omp parallel */

    return;
}


void
gradf_adjd(double *du,
           const double *x, const double *y, const double *z,
           const double *h, const size_t *sz)
{
    size_t i, j, k;
    size_t l;

    const size_t nx = sz[0];
    const size_t ny = sz[1];
    const size_t nz = sz[2];
    const size_t nxny = nx*ny;
    const size_t nxnynz = nx*ny*nz;

    const double hx = -1.0/h[0];
    const double hy = -1.0/h[1];
    const double hz = -1.0/h[2];

    /* i = 0, j = 0, k = 0 */
    l = 0;
    du[l] =
        hx*(x[l+1]-x[l]) +
        hy*(y[l+nx]-y[l]) +
        hz*(z[l+nxny]-z[l]);

#pragma omp parallel private(i,j,k,l) if (nxny*nz > 16*16*16)
{
    /* i = 0, k = 0 */
    #pragma omp for schedule(static)
    for(l = nx; l < nxny; l += nx) {
        du[l] =
            hx*(x[l+1]-x[l]) +
            hy*(y[l]-y[l-nx]) +
            hz*(z[l+nxny]-z[l]);
    }

    /* j = 0, k = 0 */
    #pragma omp for schedule(static)
    for(l = 1; l < nx; ++l) {
        du[l] =
            hx*(x[l]-x[l-1]) +
            hy*(y[l+nx]-y[l]) +
            hz*(z[l+nxny]-z[l]);
    }

    /* k = 0 */
    #pragma omp for schedule(static) collapse(2)
    for(j = nx; j < nxny; j += nx) {
        for(i = 1; i < nx; ++i) {
            l = i + j;
            du[l] =
                hx*(x[l]-x[l-1]) +
                hy*(y[l]-y[l-nx]) +
                hz*(z[l+nxny]-z[l]);
        }
    }

    /* interior loop */
    #pragma omp for schedule(static)
    for(k = nxny; k < nxnynz; k += nxny) {
        /* i = 0, j = 0 */
        l = k;
        du[l] =
            hx*(x[l+1]-x[l]) +
            hy*(y[l+nx]-y[l]) +
            hz*(z[l]-z[l-nxny]);

        /* j = 0 */
        l = 1 + k;
        for(i = 1; i < nx; ++i, ++l) {
            du[l] =
                hx*(x[l]-x[l-1]) +
                hy*(y[l+nx]-y[l]) +
                hz*(z[l]-z[l-nxny]);
        }

        for(j = nx; j < nxny; j += nx) {
            /* i = 0 */
            l = j + k;
            du[l] =
                hx*(x[l+1]-x[l]) +
                hy*(y[l]-y[l-nx]) +
                hz*(z[l]-z[l-nxny]);

            l = 1 + j + k;
            for(i = 1; i < nx; ++i, ++l) {
                du[l] =
                    hx*(x[l]-x[l-1]) +
                    hy*(y[l]-y[l-nx]) +
                    hz*(z[l]-z[l-nxny]);
            }
        }
    }

} /* omp parallel */

    return;
}
