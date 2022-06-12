#include <inttypes.h>
#include <omp.h>
#include <string.h>
#include "mx_util.h"
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

//#include <blas.h>
// #include <inttypes.h>
//#include "mex.h"
//#include "mx_blas.h"
//

extern double dnrm2_(const ptrdiff_t*, double*, const ptrdiff_t*);
extern double ddot_(const ptrdiff_t*, double*, const ptrdiff_t*, double*, const ptrdiff_t*);

extern double daxpy_(const ptrdiff_t*, double*, double*, const ptrdiff_t*, double*, const ptrdiff_t*);

extern double dscal_(const ptrdiff_t*, double*, double*, const ptrdiff_t*);




///* Level 1 BLAS */
//void
//mx_swap(mxArray *mxx, mxArray *mxy)
//{
//    const ptrdiff_t incx = 1;
//    const ptrdiff_t incy = 1;
//    const ptrdiff_t n = (const ptrdiff_t)mxGetNumberOfElements(mxx);
//
//    if (mxIsSingle(mxx)) {
//        float *x = (float *)mxGetData(mxx);
//        float *y = (float *)mxGetData(mxy);
//        sswap(&n, x, &incx, y, &incy);
//
//    } else if (mxIsDouble(mxx)) {
//        double *x = (double *)mxGetData(mxx);
//        double *y = (double *)mxGetData(mxy);
//        dswap(&n, x, &incx, y, &incy);
//
//    } else {
//        mexErrMsgTxt("not implemented");
//    }
//
//    return;
//} /* mx_swap() */
//


void
mx_scal(const double a, double *mxx, const size_t* szz)
{
    const ptrdiff_t incx = 1;
    const size_t num_szz = szz[0] * szz[1] * szz[2];
    const ptrdiff_t n = (const ptrdiff_t)num_szz;

    double *x = (mxx);

    dscal_(&n, &a, x, &incx);
    // DSCAL(&n, &a, x, &incx);

    return;
} /* mx_scal() */

//
//void
//mx_copy(const mxArray *mxx, mxArray *mxy)
//{
//    const ptrdiff_t incx = 1;
//    const ptrdiff_t incy = 1;
//    const ptrdiff_t n = (const ptrdiff_t)mxGetNumberOfElements(mxx);
//
//    if (mxIsSingle(mxx)) {
//        float *y = (float *)mxGetData(mxy);
//        const float *x = (const float *)mxGetData(mxx);
//        scopy(&n, x, &incx, y, &incy);
//
//    } else if (mxIsDouble(mxx)) {
//        double *y = (double *)mxGetData(mxy);
//        const double *x = (const double *)mxGetData(mxx);
//        dcopy(&n, x, &incx, y, &incy);
//
//    } else {
//        mexErrMsgTxt("not implemented");
//    }
//
//    return;
//} /* mx_copy() */
//


void
mx_axpy(const double a, const double* mxx, double* mxy, const size_t* szz)
{
    const ptrdiff_t incx = 1;
    const ptrdiff_t incy = 1;

    const size_t num_szz = szz[0] * szz[1] * szz[2];
    const ptrdiff_t n = num_szz;

    double* y = (mxy);
    const double* x = (mxx);
    daxpy_(&n, &a, x, &incx, y, &incy);
    // DAXPY(&n, &a, x, &incx, y, &incy);

    return;
} /* mx_axpy() */


double
mx_dot(const double*mxx, const double*mxy, const size_t* szz)
{
    double p = -1.0;

    const ptrdiff_t incx = 1;
    const ptrdiff_t incy = 1;

    const size_t num_szz = szz[0] * szz[1] * szz[2];
    const ptrdiff_t n = num_szz;

    const double* x = (mxx);
    const double* y = (mxy);

    p = ddot_(&n, x, &incx, y, &incy);
    //  p = DDOT(&n, x, &incx, y, &incy);

    return p;
} /* mx_dot() */

//


double
mx_nrm2(const double *mxx, const size_t *szz)
{
    double p = -1.0;

    const ptrdiff_t incx = 1;

    const size_t num_szz = szz[0] * szz[1] * szz[2];
    const ptrdiff_t n = (const ptrdiff_t)num_szz;

    p = dnrm2_(&n, mxx, &incx);
    // p = DNRM2(&n, mxx, &incx);

    return p;
} /* mx_nrm2() */
//
//
//double
//mx_asum(const mxArray *mxx)
//{
//    double p = -1.0;
//
//    const ptrdiff_t incx = 1;
//    const ptrdiff_t n = (const ptrdiff_t)mxGetNumberOfElements(mxx);
//
//    if (mxIsSingle(mxx)) {
//        const float *x = (const float *)mxGetData(mxx);
//        p = (double)sasum(&n, x, &incx);
//
//    } else if (mxIsDouble(mxx)) {
//        const double *x = (const double *)mxGetData(mxx);
//        p = dasum(&n, x, &incx);
//
//    } else {
//        mexErrMsgTxt("not implemented");
//    }
//
//    return p;
//} /* mx_asum() */
