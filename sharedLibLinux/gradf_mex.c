#include <inttypes.h>
#include <omp.h>
#include <stddef.h>

void
gradff(float *dx, float *dy, float *dz,
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

#pragma omp parallel private(i,j,k,l) if (nxny*nz > 16*16*16)
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
            dx[l] = hx*(u[l]-u[l-1]);
            dy[l] = hy*(u[l+nx]-u[l]);
            dz[l] = hz*(u[l+nxny]-u[l]);

        }

        /* j = ny-1 */
        l = NY + k;
        for(i = 0; i < NX; ++i, ++l) {
            dx[l] = hx*(u[l+1]-u[l]);
            dy[l] = hy*(u[l]-u[l-nx]);
            dz[l] = hz*(u[l+nxny]-u[l]);
        }

        /* i = nx-1, j = ny-1 */
        l = NX + NY + k;
        dx[l] = hx*(u[l]-u[l-1]);
        dy[l] = hy*(u[l]-u[l-nx]);
        dz[l] = hz*(u[l+nxny]-u[l]);

    }

    /* k = nz-1 */
    #pragma omp for schedule(static) collapse(2)
    for(j = 0; j < NY; j += nx) {
        for(i = 0; i < NX; ++i) {
            l = i + j + NZ;
            dx[l] = hx*(u[l+1]-u[l]);
            dy[l] = hy*(u[l+nx]-u[l]);
            dz[l] = hz*(u[l]-u[l-nxny]);
        }
    }

    /* j = ny-1, k = nz-1 */
    l = NY + NZ;
    #pragma omp for schedule(static)
    for(i = 0; i < NX; ++i) {
        l = i + NY + NZ;
        dx[l] = hx*(u[l+1]-u[l]);
        dy[l] = hy*(u[l]-u[l-nx]);
        dz[l] = hz*(u[l]-u[l-nxny]);
    }

    /* i = nx-1, k = nz-1 */
    l = NX + NZ;
    #pragma omp for schedule(static)
    for(j = 0; j < NY; j += nx) {
        l = NX + j + NZ;
        dx[l] = hx*(u[l]-u[l-1]);
        dy[l] = hy*(u[l+nx]-u[l]);
        dz[l] = hz*(u[l]-u[l-nxny]);
    }

} /* omp parallel */

    /* i = nx-1, j = ny-1, k = nz-1 */
    l = NX + NY + NZ;
    dx[l] = dx[l-1];
    dy[l] = dy[l-nx];
    dz[l] = dz[l-nxny];

    return;
}

// static int print_done = 0;

void
gradfd(double *dx, double *dy, double *dz,
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

    // if (print_done == 0) {
    //     int nz_cnt = 0;
    //     printf("JL : gradfp_adjd : Size = [%d %d %d]  H = [%.2f 0.2f 0.2f]\n", nx, ny, nz, hx, hy, hz);
    //     for (int i1 = 0; i1 < nxnynz && nz_cnt < 5; i1++) {
    //         if (x[i1] > 0) {
    //             printf("Non Zero X = %e\n", x[i1]);
    //             nz_cnt++;
    //         }
    //     }
    //     print_done = 1;
    // }

#pragma omp parallel private(i,j,k,l) if (nxny*nz > 16*16*16)
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
            dx[l] = hx*(u[l]-u[l-1]);
            dy[l] = hy*(u[l+nx]-u[l]);
            dz[l] = hz*(u[l+nxny]-u[l]);

        }

        /* j = ny-1 */
        l = NY + k;
        for(i = 0; i < NX; ++i, ++l) {
            dx[l] = hx*(u[l+1]-u[l]);
            dy[l] = hy*(u[l]-u[l-nx]);
            dz[l] = hz*(u[l+nxny]-u[l]);
        }

        /* i = nx-1, j = ny-1 */
        l = NX + NY + k;
        dx[l] = hx*(u[l]-u[l-1]);
        dy[l] = hy*(u[l]-u[l-nx]);
        dz[l] = hz*(u[l+nxny]-u[l]);

    }

    /* k = nz-1 */
    #pragma omp for schedule(static) collapse(2)
    for(j = 0; j < NY; j += nx) {
        for(i = 0; i < NX; ++i) {
            l = i + j + NZ;
            dx[l] = hx*(u[l+1]-u[l]);
            dy[l] = hy*(u[l+nx]-u[l]);
            dz[l] = hz*(u[l]-u[l-nxny]);
        }
    }

    /* j = ny-1, k = nz-1 */
    l = NY + NZ;
    #pragma omp for schedule(static)
    for(i = 0; i < NX; ++i) {
        l = i + NY + NZ;
        dx[l] = hx*(u[l+1]-u[l]);
        dy[l] = hy*(u[l]-u[l-nx]);
        dz[l] = hz*(u[l]-u[l-nxny]);
    }

    /* i = nx-1, k = nz-1 */
    l = NX + NZ;
    #pragma omp for schedule(static)
    for(j = 0; j < NY; j += nx) {
        l = NX + j + NZ;
        dx[l] = hx*(u[l]-u[l-1]);
        dy[l] = hy*(u[l+nx]-u[l]);
        dz[l] = hz*(u[l]-u[l-nxny]);
    }

} /* omp parallel */

    /* i = nx-1, j = ny-1, k = nz-1 */
    l = NX + NY + NZ;
    dx[l] = dx[l-1];
    dy[l] = dy[l-nx];
    dz[l] = dz[l-nxny];

    return;
}
