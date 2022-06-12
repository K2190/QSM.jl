#pragma once


/* Level 1 BLAS */
//extern void mx_swap(mxArray *mxx, mxArray *mxy);
extern void mx_scal(const double a, double *mxx, const size_t* szz);
//extern void mx_copy(const mxArray *mxx, mxArray *mxy);
extern void mx_axpy(const double a, const double *mxx, double *mxy, const size_t* szz);
//
extern double mx_dot(const double*mxx, const double*mxy, const size_t* szz);
extern double mx_nrm2(const double *mxx, const size_t* szz);
//extern double mx_asum(const mxArray *mxx);
