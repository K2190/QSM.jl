#include <inttypes.h>
#include <omp.h>
#include <stdlib.h>
#include "boundary_mask_mex.h"
#include <stdio.h>

void
boundary_mask(uint8_t *B, const uint8_t *G, const size_t *sz)
{
    size_t i, j, k;
    size_t l;

    const size_t nx = sz[0];
    const size_t ny = sz[1];
    const size_t nz = sz[2];
    const size_t nxny = nx*ny;

    const size_t NX = nx-2;
    const size_t NY = nx*(ny-2);
    const size_t NZ = nxny*(nz-2);

    const size_t nx2 = 2*nx;
    const size_t nxny2 = 2*nxny;
    
    // printf("INSIDE boundary_mask: nx= %d  ny= %d  nz= %d\n", nx, ny, nz);

    /* offset indices */
    const size_t o110 = 1 + nx +    0;
    const size_t o101 = 1 +  0 + nxny;
    const size_t o011 = 0 + nx + nxny;
    const size_t o111 = 1 + nx + nxny;

    uint8_t *b = (uint8_t *)calloc(nx*ny*nz, sizeof(*G));

    /* boundary of grid */
    #pragma omp parallel for private(i,j,k,l) schedule(static) \
        if (nxny*nz > 32*32*32)
    for(k = nxny; k <= NZ; k += nxny) {
        for(j = nx; j <= NY; j += nx) {
            l = 1 + j + k;
            for(i = 1; i <= NX; ++i, ++l) {
                if ((i == 1) || (j == nx) || (k == nxny) ||
                        (i == NX) || (j == NY) || (k == NZ)) {
                    b[l] = G[l];
                }
            }
        }
    }

    /* interior */
    #pragma omp parallel for private(i,j,k,l) schedule(static) \
        if (nxny*nz > 32*32*32)
    for(k = nxny; k <= NZ; k += nxny2) {
        for(j = nx; j <= NY; j += nx2) {
            l = 1 + j + k;
            for(i = 1; i <= NX; i += 2, l += 2) {

                if (!(G[l] && G[l+1] && G[l+nx] && G[l+nxny] &&
                        G[l+o110] && G[l+o101] && G[l+o011] && G[l+o111])) {

                    b[l] = G[l];
                    b[l+1] = G[l+1];
                    b[l+nx] = G[l+nx];
                    b[l+nxny] = G[l+nxny];
                    b[l+o110] = G[l+o110];
                    b[l+o101] = G[l+o101];
                    b[l+o011] = G[l+o011];
                    b[l+o111] = G[l+o111];
                }
            }
        }
    }

    /* grow boundary band */
    #pragma omp parallel for private(i,j,k,l) schedule(static) \
        if (nxny*nz > 32*32*32)
    for(k = nxny; k <= NZ; k += nxny) {
        for(j = nx; j <= NY; j += nx) {
            l = 1 + j + k;
            for(i = 1; i <= NX; ++i, ++l) {
                if (G[l]) {
                    B[l] = b[l-nxny] || b[l-nx] || b[l-1] || b[l] ||
                        b[l+1] || b[l+nx] || b[l+nxny];
                }
            }
        }
    }

    if (NULL != b) {
        free(b);
        b = NULL;
    }

    return;
}
