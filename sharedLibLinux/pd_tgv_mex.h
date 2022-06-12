#pragma once

#include <inttypes.h>
#include <stddef.h>

struct PD_TGV
{
    /* primal */
    void* chi;
    void* eta;
    void* w;
    void* chi_;
    void* eta_;
    void* w_;
    /* dual */
    void* nu;
    void* p;
    void* q;
    /* misc */
    void* f;
    void* m;
    void* mi;
    void* mb;       /* unused */
    void* chi0;
    void* eta0;
    /* TODO: w0 */
};

extern void update_duald(struct PD_TGV* pd,
    const double alpha0, const double alpha1, const double sigma,
    const double* h, const size_t* sz);

extern void update_dualf(struct PD_TGV* pd,
    const double alpha0, const double alpha1, const double sigma,
    const double* h, const size_t* sz);

extern void update_primald(struct PD_TGV* pd,
    const double tau, const double* h, const size_t* sz);

extern void update_primalf(struct PD_TGV* pd,
    const double tau, const double* h, const size_t* sz);

extern uint8_t convergence_checkd(double* nr1, double* nr2,
    const double* chi, const double* chi0,
    const double* eta, const double* eta0,
    const double tol, const size_t N);

extern uint8_t convergence_checkf(double* nr1, double* nr2,
    const float* chi, const float* chi0,
    const float* eta, const float* eta0,
    const double tol, const size_t N);

