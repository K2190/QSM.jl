#include <inttypes.h>
#include <omp.h>
#include <stddef.h>
#include "coarsen_grid_mex.h"
#include <stdio.h>

void
coarsen_grid(uint8_t *G2, const uint8_t *G, const size_t *sz2, const size_t *sz)
{
    size_t i2, j2, k2;
    size_t l2, lk2;
    size_t i, j, k;
    size_t l, lk;

    uint8_t bi, bj, bk;

    const size_t nx = sz[0];
    const size_t ny = sz[1];
    const size_t nz = sz[2];
    const size_t nxny = nx*ny;

    const size_t nx2 = sz2[0];
    const size_t ny2 = sz2[1];
    const size_t nz2 = sz2[2];
    const size_t nxny2 = nx2*ny2;

    const size_t NX = nx-2;
    const size_t NY = ny-2;
    const size_t NZ = nz-2;

    const size_t NX2 = nx2-1;
    const size_t NY2 = ny2-1;
    const size_t NZ2 = nz2-1;

    // printf("INSIDE coarsen_grid: nx= %d  ny= %d  nz= %d\n", nx, ny, nz);
    // printf("INSIDE coarsen_grid: nx2= %d  ny2= %d  nz2= %d\n", nx2, ny2, nz2);

    /* offset indices */
    const size_t o110 = 1 + nx +    0;
    const size_t o101 = 1 +  0 + nxny;
    const size_t o011 = 0 + nx + nxny;
    const size_t o111 = 1 + nx + nxny;


    #pragma omp parallel for private(l2) schedule(static) \
        if (nxny2*nz2 > 64*64*64)
    for(l2 = 0; l2 < nxny2*nz2; ++l2) {
        G2[l2] = 0;
    }

    #pragma omp parallel for private(i2,j2,k2,lk2,l2,i,j,k,lk,l,bi,bj,bk) \
        schedule(static)
    for(k2 = 1; k2 < NZ2; ++k2) {
        k = (k2<<1)-1;
        lk = nxny*k;
        lk2 = nxny2*k2;

        for(j2 = 1; j2 < NY2; ++j2) {
            j = (j2<<1)-1;
            l = 1 + nx*j + lk;
            l2 = 1 + nx2*j2 + lk2;

            for(i2 = 1; i2 < NX2; ++i2, ++l2, l += 2) {
                i = (i2<<1)-1;

                bi = !(i < NX);
                bj = !(j < NY);
                bk = !(k < NZ);

                G2[l2] = G[l] &&
                    (bi || G[l+1]) &&
                    (bj || G[l+nx]) &&
                    (bk || G[l+nxny]) &&
                    (bi || bj || G[l+o110]) &&
                    (bi || bk || G[l+o101]) &&
                    (bj || bk || G[l+o011]) &&
                    (bi || bj || bk || G[l+o111]);
            }
        }
    }

    return;
}
