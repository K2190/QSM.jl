#include <inttypes.h>
#include <memory.h>
#include "boundary_mask_l1_mex.h"
#include <stdio.h>
#include <stdlib.h>

void
boundary_mask_l1(uint8_t *B,
                 const uint8_t *G, const size_t *sz, const uint8_t l1)
{
    size_t i, j, k;
    size_t l;

    const size_t nx = sz[0];
    const size_t ny = sz[1];
    const size_t nz = sz[2];
    const size_t nxny = nx*ny;
    const size_t nxnynz = nx*ny*nz;

    size_t NX = nx-1;
    size_t NY = nx*(ny-1);
    size_t NZ = nxny*(nz-1);

    uint8_t *GL = (uint8_t *)calloc(nx*ny*nz, sizeof(*G));

    printf("Inside boundary_mask_l1...\n");
    printf("Size = %d  %d  %d\n", sz[0], sz[1], sz[2]);
    printf("l1 = %d\n", l1);

    if (GL == NULL) {
        printf("NO MEMORY AVAILABLE..........\n");
        return;
    }


    #pragma omp parallel for private(i,j,k,l) schedule(static) \
        if(nxny*nz > 32*32*32)
    for(k = 0; k < nxnynz; k += nxny) {
        for(j = 0; j < nxny; j += nx) {
            l = j + k;
            for(i = 0; i < nx; ++i, ++l) {

                if ((i == 0) || (j == 0) || (k == 0)
                        || (i == NX) || (j == NY) || (k == NZ)) {

                    B[l] = GL[l] = G[l];

                } else {
                    B[l] = GL[l] = G[l] &&
                        !(G[l-1] && G[l+1]
                        && G[l-nx] && G[l+nx]
                        && G[l-nxny] && G[l+nxny]);
                }
            }
        }
    }

    for(size_t L1 = 1; L1 < l1; ++L1) {
        #pragma omp parallel for private(i,j,k,l) schedule(static) \
            if(nxny*nz > 32*32*32)
        for(k = L1*nxny; k < nxnynz-L1*nxny; k += nxny) {
            for(j = L1*nx; j < nxny-L1*nx; j += nx) {
                l = L1 + j + k;
                for(i = L1; i < nx-L1; ++i, ++l) {

                    if (G[l] && !GL[l] &&
                        ((GL[l-1] == L1) || (GL[l+1] == L1) ||
                        (GL[l-nx] == L1) || (GL[l+nx] == L1) ||
                        (GL[l-nxny] == L1) || (GL[l+nxny] == L1))) {

                        B[l] = 1;
                        GL[l] = L1+1;
                    }
                }
            }
        }
    }

    if (NULL != GL) {
        free(GL);
        GL = NULL;
    }

    return;
}
