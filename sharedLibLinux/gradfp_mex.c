#include <inttypes.h>
#include <omp.h>
#include <stddef.h>

static int prdone = 0;

void
gradfpf(float *dx, float *dy, float *dz,
        const float *u, const double *h, const size_t *sz)
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

    const float hx = (float)(1.0/h[0]);
    const float hy = (float)(1.0/h[1]);
    const float hz = (float)(1.0/h[2]);

    if(prdone < 3) {
        printf("INSIDE gradfpf\n");
        int cnt = 0;

        printf("Size = [%d %d %d]   H = [%f %f %f]\n", nx, ny, nz, hx, hy, hz);
        const size_t nxnynz = nx*ny*nz;
        for(size_t i1 = 0; i1 < nxnynz && cnt < 10; i1++) {
            if (u[i1] > 0) {
                printf("Data = %e\n", u[i1]);
                cnt++;
            }
        }
        prdone++;
    }

#pragma omp parallel private(i,j,k,l)
{
    #pragma omp for schedule(static)
    for(k = 0; k < NZ; k += nxny) {
        for(j = 0; j < NY; j += nx) {
            l = j + k;
            for(i = 0; i < NX; ++i, ++l) {
                dx[l] = hx*(u[l+1]-u[l]);
                dy[l] = hy*(u[l+nx]-u[l]);
                dz[l] = hz*(u[l+nxny]-u[l]);
            }

            /* i = nx-1 */
            l = NX + j + k;
            dx[l] = hx*(u[l-NX]-u[l]);
            dy[l] = hy*(u[l+nx]-u[l]);
            dz[l] = hz*(u[l+nxny]-u[l]);

        }

        /* j = ny-1 */
        l = NY + k;
        for(i = 0; i < NX; ++i, ++l) {
            dy[l] = hy*(u[l-NY]-u[l]);
            dx[l] = hx*(u[l+1]-u[l]);
            dz[l] = hz*(u[l+nxny]-u[l]);
        }

        /* i = nx-1, j = ny-1 */
        l = NX + NY + k;
        dy[l] = hy*(u[l-NY]-u[l]);
        dx[l] = hx*(u[l-NX]-u[l]);
        dz[l] = hz*(u[l+nxny]-u[l]);

    }

    /* k = nz-1 */
    #pragma omp for schedule(static) collapse(2)
    for(j = 0; j < NY; j += nx) {
        for(i = 0; i < NX; ++i) {
            l = i + j + NZ;
            dz[l] = hz*(u[l-NZ]-u[l]);
            dx[l] = hx*(u[l+1]-u[l]);
            dy[l] = hy*(u[l+nx]-u[l]);
        }
    }

    /* j = ny-1, k = nz-1 */
    #pragma omp for schedule(static)
    for(i = 0; i < NX; ++i) {
        l = i + NY + NZ;
        dz[l] = hz*(u[l-NZ]-u[l]);
        dy[l] = hy*(u[l-NY]-u[l]);
        dx[l] = hx*(u[l+1]-u[l]);
    }

    /* i = nx-1, k = nz-1 */
    #pragma omp for schedule(static)
    for(j = 0; j < NY; j += nx) {
        l = j + NX + NZ;
        dz[l] = hz*(u[l-NZ]-u[l]);
        dx[l] = hx*(u[l-NX]-u[l]);
        dy[l] = hy*(u[l+nx]-u[l]);
    }

} /* omp parallel */

    /* i = nx-1, j = ny-1, k = nz-1 */
    l = NX + NY + NZ;
    dz[l] = hz*(u[l-NZ]-u[l]);
    dy[l] = hy*(u[l-NY]-u[l]);
    dx[l] = hx*(u[l-NX]-u[l]);

    return;
}

void
gradfpd(double *dx, double *dy, double *dz,
        const double *u, const double *h, const size_t *sz)
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

    const double hx = 1.0/h[0];
    const double hy = 1.0/h[1];
    const double hz = 1.0/h[2];

    
    // if(prdone < 3) {
    //     printf("INSIDE gradfpd\n");
    //     int cnt = 0;

    //     printf("Size = [%d %d %d]   H = [%f %f %f]\n", nx, ny, nz, hx, hy, hz);
    //     const size_t nxnynz = nx*ny*nz;
    //     for(size_t i1 = 0; i1 < nxnynz && cnt < 10; i1++) {
    //         if (u[i1] > 0) {
    //             printf("Data = %e\n", u[i1]);
    //             cnt++;
    //         }
    //     }
    //     printf("PRDONE ==== %d\n", prdone);
    //     prdone++;
    // }


#pragma omp parallel private(i,j,k,l)
{
    #pragma omp for schedule(static)
    for(k = 0; k < NZ; k += nxny) {
        for(j = 0; j < NY; j += nx) {
            l = j + k;
            for(i = 0; i < NX; ++i, ++l) {
                dx[l] = hx*(u[l+1]-u[l]);
                dy[l] = hy*(u[l+nx]-u[l]);
                dz[l] = hz*(u[l+nxny]-u[l]);
            }

            /* i = nx-1 */
            l = NX + j + k;
            dx[l] = hx*(u[l-NX]-u[l]);
            dy[l] = hy*(u[l+nx]-u[l]);
            dz[l] = hz*(u[l+nxny]-u[l]);

        }

        /* j = ny-1 */
        l = NY + k;
        for(i = 0; i < NX; ++i, ++l) {
            dy[l] = hy*(u[l-NY]-u[l]);
            dx[l] = hx*(u[l+1]-u[l]);
            dz[l] = hz*(u[l+nxny]-u[l]);
        }

        /* i = nx-1, j = ny-1 */
        l = NX + NY + k;
        dy[l] = hy*(u[l-NY]-u[l]);
        dx[l] = hx*(u[l-NX]-u[l]);
        dz[l] = hz*(u[l+nxny]-u[l]);

    }

    /* k = nz-1 */
    #pragma omp for schedule(static) collapse(2)
    for(j = 0; j < NY; j += nx) {
        for(i = 0; i < NX; ++i) {
            l = i + j + NZ;
            dz[l] = hz*(u[l-NZ]-u[l]);
            dx[l] = hx*(u[l+1]-u[l]);
            dy[l] = hy*(u[l+nx]-u[l]);
        }
    }

    /* j = ny-1, k = nz-1 */
    #pragma omp for schedule(static)
    for(i = 0; i < NX; ++i) {
        l = i + NY + NZ;
        dz[l] = hz*(u[l-NZ]-u[l]);
        dy[l] = hy*(u[l-NY]-u[l]);
        dx[l] = hx*(u[l+1]-u[l]);
    }

    /* i = nx-1, k = nz-1 */
    #pragma omp for schedule(static)
    for(j = 0; j < NY; j += nx) {
        l = j + NX + NZ;
        dz[l] = hz*(u[l-NZ]-u[l]);
        dx[l] = hx*(u[l-NX]-u[l]);
        dy[l] = hy*(u[l+nx]-u[l]);
    }

} /* omp parallel */

    /* i = nx-1, j = ny-1, k = nz-1 */
    l = NX + NY + NZ;
    dz[l] = hz*(u[l-NZ]-u[l]);
    dy[l] = hy*(u[l-NY]-u[l]);
    dx[l] = hx*(u[l-NX]-u[l]);

    // if (prdone == 2) {
    //     int cnt = 0;

    //     const size_t nxnynz = nx*ny*nz;
    //     for(size_t i1 = 0; i1 < nxnynz && cnt < 20; i1++) {
    //         if (dx[i1] > 0) {
    //             printf("dx = %e\n", dx[i1]);
    //             cnt++;
    //         }
    //     }
    // }

    return;
}
