#include <inttypes.h>
#include <omp.h>
#include <string.h>
#include "mx_util.h"
#include <stdio.h>
#include <stdlib.h>

uint8_t*
mx_pad_boundary_uint(const size_t* sz, const uint8_t* mxx, size_t* new_sz)
{
    const size_t szp[3] = { sz[0] + 2, sz[1] + 2, sz[2] + 2 };

    new_sz[0] = szp[0];
    new_sz[1] = szp[1];
    new_sz[2] = szp[2];

    const size_t nx = sz[0];
    const size_t ny = sz[1];
    const size_t nz = sz[2];
    const size_t nxny = nx * ny;

    const size_t nxp = szp[0];
    const size_t nyp = szp[1];
    const size_t nxnyp = nxp * nyp;

    uint8_t* x = mxx;
    uint8_t* xp = (uint8_t*)malloc(szp[0] * szp[1] * szp[2] * sizeof(uint8_t));
    memset(xp, 0, szp[0] * szp[1] * szp[2]);


#pragma omp parallel for schedule(static) collapse(3) \
            if(nxny*nz > 32*32*32)
    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                xp[(i + 1) + nxp * (j + 1) + nxnyp * (k + 1)] = x[i + nx * j + nxny * k];
            }
        }
    }

    return xp;
}


double*
mx_pad_boundary(const size_t* sz, const double* mxx, size_t* new_sz)
{
    const size_t szp[3] = { sz[0] + 2, sz[1] + 2, sz[2] + 2 };

    new_sz[0] = szp[0];
    new_sz[1] = szp[1];
    new_sz[2] = szp[2];

    const size_t nx = sz[0];
    const size_t ny = sz[1];
    const size_t nz = sz[2];
    const size_t nxny = nx * ny;

    const size_t nxp = szp[0];
    const size_t nyp = szp[1];
    const size_t nxnyp = nxp * nyp;

    double* x = mxx;
    double* xp = (double*)malloc(szp[0] * szp[1] * szp[2] * sizeof(double));
    memset(xp, 0, szp[0] * szp[1] * szp[2]);


#pragma omp parallel for schedule(static) collapse(3) \
            if(nxny*nz > 32*32*32)
    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                xp[(i + 1) + nxp * (j + 1) + nxnyp * (k + 1)] = x[i + nx * j + nxny * k];
            }
        }
    }

    return xp;
}


double*
mx_unpad_boundary(const double* mxxp, size_t* mxxp_sz, double* mxx)
{

    //const size_t *szp = (const size_t *)mxGetDimensions(mxxp);
    const size_t szp[3] = { mxxp_sz[0], mxxp_sz[1], mxxp_sz[2]};
    const size_t sz[3] = { mxxp_sz[0]-2, mxxp_sz[1]-2, mxxp_sz[2]-2};

    // mxCreateNumericArray(3, sz, mxGetClassID(mxxp), mxREAL);
    
    // Get mxx as argument
    //double* mxx = (double*)malloc(sz[0] * sz[1] * sz[2] * sizeof(double));

    const size_t nx = sz[0];
    const size_t ny = sz[1];
    const size_t nz = sz[2];
    const size_t nxny = nx*ny;

    const size_t nxp = szp[0];
    const size_t nyp = szp[1];
    const size_t nxnyp = nxp*nyp;

    //if (mxIsSingle(mxxp)) {
    //    float *x = (float *)mxGetData(mxx);
    //    float *xp = (float *)mxGetData(mxxp);

    //    #pragma omp parallel for schedule(static) collapse(3) \
    //        if(nxny*nz > 32*32*32)
    //    for (size_t k = 0; k < nz; ++k) {
    //        for (size_t j = 0; j < ny; ++j) {
    //            for (size_t i = 0; i < nx; ++i) {
    //                x[i + nx*j + nxny*k] = xp[(i+1) + nxp*(j+1) + nxnyp*(k+1)];
    //            }
    //        }
    //    }

    //} else if (mxIsDouble(mxxp)) {
        double *x = (double *)mxx;
        double *xp = (double *)mxxp;

        #pragma omp parallel for schedule(static) collapse(3) \
            if(nxny*nz > 32*32*32)
        for (size_t k = 0; k < nz; ++k) {
            for (size_t j = 0; j < ny; ++j) {
                for (size_t i = 0; i < nx; ++i) {
                    x[i + nx*j + nxny*k] = xp[(i+1) + nxp*(j+1) + nxnyp*(k+1)];
                }
            }
        }

    //} else if (mxIsLogical(mxxp)) {
    //    uint8_t *x = (uint8_t *)mxGetData(mxx);
    //    uint8_t *xp = (uint8_t *)mxGetData(mxxp);

    //    #pragma omp parallel for schedule(static) collapse(3) \
    //        if(nxny*nz > 32*32*32)
    //    for (size_t k = 0; k < nz; ++k) {
    //        for (size_t j = 0; j < ny; ++j) {
    //            for (size_t i = 0; i < nx; ++i) {
    //                x[i + nx*j + nxny*k] = xp[(i+1) + nxp*(j+1) + nxnyp*(k+1)];
    //            }
    //        }
    //    }
    //}

    return mxx;
}


void
mx_zero(double* mxx)
{
    //const size_t n = mxGetNumberOfElements(mxx);
    //const size_t b = mxGetElementSize(mxx);
    //memset(mxGetData(mxx), 0, n*b);
    return;
}
