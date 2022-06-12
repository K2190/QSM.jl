#ifndef _MATLAB_MEX_H_
#define _MATLAB_MEX_H_

#include <inttypes.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stddef.h>

extern void lapwf(float* du, const float* u, const double* h, const size_t* sz);
extern void lapwd(double* du, const double* u, const double* h, const size_t* sz);
extern void gradfp_adjf(float* du,
    const float* x, const float* y, const float* z,
    const double* h, const size_t* sz);
extern void gradfp_adjd(double* du,
    const double* x, const double* y, const double* z,
    const double* h, const size_t* sz);

extern void gradfpf(float* dx, float* dy, float* dz,
    const float* u, const double* h, const size_t* sz);

extern void gradfpd(double* dx, double* dy, double* dz,
    const double* u, const double* h, const size_t* sz);



// utils/fd

extern void gradb_adjf(float* du,
    const float* x, const float* y, const float* z,
    const double* h, const size_t* sz);

extern void gradb_adjd(double* du,
    const double* x, const double* y, const double* z,
    const double* h, const size_t* sz);

extern void gradbf(float* dx, float* dy, float* dz,
    const float* u, const double* h, const size_t* sz);

extern void gradbd(double* dx, double* dy, double* dz,
    const double* u, const double* h, const size_t* sz);

extern void gradbm_adjf(float* du,
    const float* x, const float* y, const float* z,
    const uint8_t* G, const double* h, const size_t* sz);

extern void gradbm_adjd(double* du,
    const double* x, const double* y, const double* z,
    const uint8_t* G, const double* h, const size_t* sz);

extern void gradcf(float* dx, float* dy, float* dz,
    const float* u, const double* h, const size_t* sz);

extern void gradcd(double* dx, double* dy, double* dz,
    const double* u, const double* h, const size_t* sz);

extern void gradcmf(float* dx, float* dy, float* dz,
    const float* u, const uint8_t* G,
    const double* h, const size_t* sz);

extern void gradcmd(double* dx, double* dy, double* dz,
    const double* u, const uint8_t* G,
    const double* h, const size_t* sz);

extern void gradf_adjf(float* du,
    const float* x, const float* y, const float* z,
    const double* h, const size_t* sz);

extern void gradf_adjd(double* du,
    const double* x, const double* y, const double* z,
    const double* h, const size_t* sz);

extern  void gradff(float* dx, float* dy, float* dz,
    const float* u, const double* h, const size_t* sz);

void gradfd(double* dx, double* dy, double* dz,
    const double* u, const double* h, const size_t* sz);

extern  void gradfm_adjf(float* du,
    const float* x, const float* y, const float* z,
    const uint8_t* G, const double* h, const size_t* sz);

extern  void gradfm_adjd(double* du,
    const double* x, const double* y, const double* z,
    const uint8_t* G, const double* h, const size_t* sz);

extern  void gradfmf(float* dx, float* dy, float* dz,
    const float* u, const uint8_t* G,
    const double* h, const size_t* sz);

extern  void gradfmd(double* dx, double* dy, double* dz,
    const double* u, const uint8_t* G,
    const double* h, const size_t* sz);

extern  void lapf(float* du,
    const float* u, const uint8_t* G,
    const double* h, const size_t* sz);

extern  void lapd(double* du,
    const double* u, const uint8_t* G,
    const double* h, const size_t* sz);

extern  void mgpcg(const uint8_t* szv, double* v,
    const uint8_t* szf, double* f,
    const uint8_t* szm, uint8_t* mask,
    const uint8_t* szz, double* vsz,
    double tol, 
    int32_t maxit,
    double mgopts_tol,
    int32_t mgopts_maxit,
    int32_t mgopts_mu,
    int32_t mgopts_npre,
    int32_t mgopts_npost,
    int32_t mgopts_nboundary,
    int32_t mgopts_nlevels,
    int32_t verbose,
    const uint8_t* szz_out,
    double* out_mat);


#endif

