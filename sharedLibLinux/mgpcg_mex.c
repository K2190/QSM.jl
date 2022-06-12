#include <inttypes.h>
#include <stddef.h>
#include "boundary_mask_mex.h"
#include "coarsen_grid_mex.h"
#include "correct_mex.h"
#include "gauss_seidel_mex.h"
#include "lapmg_mex.h"
#include "mx_blas.h"
#include "mx_util.h"
#include "prolong_mex.h"
#include "restrict_mex.h"
#include "residual_mex.h"
#include "matlab_mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define MAX_DEPTH 16

struct MG
{
    double* f;
    double* v;
    double* r;
    uint8_t* G;
    uint8_t* B;
    double* h;
    size_t szf[3];
    size_t szv[3];
    size_t szr[3];
    size_t szG[3];
};

/* TODO: return convergence history struct */
struct PCG
{
    uint8_t  isconverged;
    int32_t  iters;
    double *resnorm;
    double   reltol;
    const double tol;
    const int32_t maxit;
};

struct Iter
{
    const int32_t mg;
    const int32_t mu;
    const int32_t pre;
    const int32_t post;
    const int32_t bs;
    const int32_t l;
    const double tol;
};

static void
mx_cycle(struct MG* mg, const struct Iter* iter, const int32_t l)
{
    // printf("mx_cycle : ITER_L = %d   L = %d\n", iter->l, l);

    if (l < iter->l - 1) {

        // printf("INSIDE mx_cycle .... IF LOOP...l = %d\n", l);

        const int32_t bs = (1 << l) * iter->bs;


        double tempTol1 = mx_nrm2(mg[l].v, (const size_t*)mg[l].szv);
        printf(">>>>>> 11111 >>>>>>>INSIDE mx_cycle .... IF LOOP...tempTol1 = %f\n", tempTol1);
        tempTol1 = mx_nrm2(mg[l].f, (const size_t*)mg[l].szf);
        printf(">>>>>> 11111 >>>>>>INSIDE mx_cycle .... IF LOOP...tempTol1 = %f\n", tempTol1);


        gauss_seideld(mg[l].v, mg[l].f, mg[l].G, mg[l].h, iter->pre, 1, mg[l].szf);
        

        tempTol1 = mx_nrm2(mg[l].v, (const size_t*)mg[l].szv);
        printf(">>>>>> 22222 >>>>>>>INSIDE mx_cycle .... IF LOOP...tempTol1 = %f\n", tempTol1);
        tempTol1 = mx_nrm2(mg[l].f, (const size_t*)mg[l].szf);
        printf(">>>>>> 22222 >>>>>>INSIDE mx_cycle .... IF LOOP...tempTol1 = %f\n", tempTol1);


        gauss_seideld(mg[l].v, mg[l].f, mg[l].B, mg[l].h,  bs, 1, mg[l].szf);


        // tempTol1 = mx_nrm2(mg[l].v, (const size_t*)mg[l].szv);
        // printf(">>>>>> 33333 >>>>>>>INSIDE mx_cycle .... IF LOOP...tempTol1 = %f\n", tempTol1);
        // tempTol1 = mx_nrm2(mg[l].f, (const size_t*)mg[l].szf);
        // printf(">>>>>> 33333 >>>>>>INSIDE mx_cycle .... IF LOOP...tempTol1 = %f\n", tempTol1);


        //mx_residual(mg[l].r, mg[l].f, mg[l].v, mg[l].G, mg[l].h);
        residuald(mg[l].r, mg[l].f, mg[l].v, mg[l].G, mg[l].h, mg[l].szf);

        // tempTol1 = mx_nrm2(mg[l].v, (const size_t*)mg[l].szv);
        // printf(">>>>>> 44444 >>>>>>>INSIDE mx_cycle .... IF LOOP...tempTol1 = %f\n", tempTol1);
        // tempTol1 = mx_nrm2(mg[l].f, (const size_t*)mg[l].szf);
        // printf(">>>>>> 44444 >>>>>>INSIDE mx_cycle .... IF LOOP...tempTol1 = %f\n", tempTol1);

    // int sample_cnt = 0;
    // uint8_t* dG = mg[l + 1].G;
    // int nZZZ = mg[l + 1].szG[0] * mg[l + 1].szG[1] * mg[l + 1].szG[2];
    // for (int i = 0; i < nZZZ; i++) {
       // if (dG[i] > 0) {
           // sample_cnt++;
       // }
    // }
    // printf("Mask is non zero for mg[l + 1].G %d\n", sample_cnt);


        //mx_restrict(mg[l + 1].f, mg[l].r, mg[l + 1].G);
        restrictd(mg[l + 1].f, mg[l].r, mg[l + 1].G, mg[l + 1].szf, mg[l].szr);


        // tempTol1 = mx_nrm2(mg[l].v, (const size_t*)mg[l].szv);
        // printf(">>>>>> 55555 >>>>>>>INSIDE mx_cycle .... IF LOOP...tempTol1 = %f\n", tempTol1);
        // tempTol1 = mx_nrm2(mg[l].f, (const size_t*)mg[l].szf);
        // printf(">>>>>> 55555 >>>>>>INSIDE mx_cycle .... IF LOOP...tempTol1 = %f\n", tempTol1);


        int num_szv_1 = mg[l + 1].szv[0] * mg[l + 1].szv[1] * mg[l + 1].szv[2];
        for (int kk = 0; kk < num_szv_1; kk++)
            mg[l + 1].v[kk] = 0;

        for (int32_t i = 0; i < iter->mu; ++i) {
            // printf("INSIDE mx_cycle ..iter->mu  FOR LOOP...calling mx_cycle at I = %d\n", i);

            mx_cycle(mg, iter, l + 1);
            if (l + 1 == iter->l - 1) { 
                // printf("INSIDE mx_cycle... IF LOOP...calling break at I = %d  (l+1) = %d \n", i, l+1);
                break;
            }
        }

        correctd(mg[l].v, mg[l + 1].v, mg[l].G, mg[l].szv, mg[l + 1].szv);

        // tempTol1 = mx_nrm2(mg[l].v, (const size_t*)mg[l].szv);
        // printf(">>>>>> 66666 >>>>>>>INSIDE mx_cycle .... IF LOOP...tempTol1 = %f\n", tempTol1);
        // tempTol1 = mx_nrm2(mg[l].f, (const size_t*)mg[l].szf);
        // printf(">>>>>> 66666 >>>>>>INSIDE mx_cycle .... IF LOOP...tempTol1 = %f\n", tempTol1);


        gauss_seideld(mg[l].v, mg[l].f, mg[l].B, mg[l].h,  bs, 0, mg[l].szf);

        // tempTol1 = mx_nrm2(mg[l].v, (const size_t*)mg[l].szv);
        // printf(">>>>>> 77777 >>>>>>>INSIDE mx_cycle .... IF LOOP...tempTol1 = %f\n", tempTol1);
        // tempTol1 = mx_nrm2(mg[l].f, (const size_t*)mg[l].szf);
        // printf(">>>>>> 77777 >>>>>>INSIDE mx_cycle .... IF LOOP...tempTol1 = %f\n", tempTol1);



        gauss_seideld(mg[l].v, mg[l].f, mg[l].G, mg[l].h, iter->post, 0, mg[l].szf);


        // tempTol1 = mx_nrm2(mg[l].v, (const size_t*)mg[l].szv);
        // printf(">>>>>> 88888 >>>>>>>INSIDE mx_cycle .... IF LOOP...tempTol1 = %f\n", tempTol1);
        // tempTol1 = mx_nrm2(mg[l].f, (const size_t*)mg[l].szf);
        // printf(">>>>>> 88888 >>>>>>INSIDE mx_cycle .... IF LOOP...tempTol1 = %f\n", tempTol1);


    }
    else {

        const double tol =  1e-8;  // CHECK THIS ????????
        const double reltol = tol * mx_nrm2(mg[l].f, (const size_t*)mg[l].szf);

        // printf("mx_cycle  ELSE LOOP before FOR LOOP STARTS reltol = %f\n", reltol);

        for (int32_t i = 0; i < 1024; ++i) {

            // double tempTolQ = mx_nrm2(mg[l].v, mg[l].szv);
            // printf(">>> QQQQQ >>>>>>>>>> tempTol1 = %f   L = %d\n", tempTolQ, l);
            // tempTolQ = mx_nrm2(mg[l].f, mg[l].szf);
            // printf(">>> QQQQQ >>>>>>>>> tempTol1 = %f\n", tempTolQ);


            gauss_seideld(mg[l].v, mg[l].f, mg[l].G, mg[l].h, 128, 1, mg[l].szf);
            gauss_seideld(mg[l].v, mg[l].f, mg[l].G, mg[l].h,  128, 0, mg[l].szf);

            residuald(mg[l].r, mg[l].f, mg[l].v, mg[l].G, mg[l].h, mg[l].szf);

            // double tempTol1 = mx_nrm2(mg[l].v, mg[l].szv);
            // printf(">>> 99999 >>>>>>>>>> tempTol1 = %f\n", tempTol1);
            // tempTol1 = mx_nrm2(mg[l].f, mg[l].szf);
            // printf(">>> 99999 >>>>>>>>> tempTol1 = %f\n", tempTol1);

            if (mx_nrm2(mg[l].r, (const size_t*)mg[l].szr) <= reltol) {
                // printf("INSIDE mx_cycle .... CALLING BREAK AT I = %d\n", i);
                break;
            }
        }
    }

    return;
}



static void
mx_mg(struct MG* mg, const struct Iter* iter)
{
    if (iter->tol >= 0.0) {
        const double reltol = iter->tol * mx_nrm2(mg[0].f, mg[0].szf);

        for (int32_t i = 0; i < iter->mg; ++i) {
            residuald(mg[0].r, mg[0].f, mg[0].v, mg[0].G, mg[0].h, mg[0].szf);

            if (mx_nrm2(mg[0].r, (const size_t*)mg[0].szr) < reltol) {
                // printf("INSIDE mx_mg... FOR LOOP...calling break at I = %d \n", i);
                break;
            }

            mx_cycle(mg, iter, 0);
        }

    }
    else {
        for (int32_t i = 0; i < iter->mg; ++i) {
            // printf("INSIDE mx_mg ELSE LOOP...calling mx_cycle at I = %d\n", i);
            mx_cycle(mg, iter, 0);
        }
    }

    return;
}

static void
mx_pcg(double* x, struct PCG* ch, struct MG* mg,
    const struct Iter* iter, const uint8_t verbose,
    const size_t* szv,
    const size_t* szf,
    const size_t* szG,
    const size_t* szz)
{
    const double ONE = 1.0;

    double* resnorm = ch->resnorm;

    /* initialize pcg scalars */
    double a = 1.0;
    double beta = 1.0;
    double rho = 1.0;
    double prevrho = rho;

    double res = mx_nrm2(mg[0].f, (const size_t*)mg[0].szf);
    const double reltol = ch->tol * res;

    // printf("INSIDE mx_pcg AT BEGINING : res from mx_nrm2 = %f\n", res);
    // printf("INSIDE mx_pcg AT BEGINING: reltol = %f\n", reltol);

    // printf("F Matrix inside mgpcg_init...\n");
    // printf("f_matrix at 22149 = %f\n", mg[0].f[22149]);
    // printf("f_matrix at 22150 = %f\n", mg[0].f[22150]);
    // printf("f_matrix at 22151 = %f\n", mg[0].f[22151]);
    // printf("f_matrix at 22152 = %f\n", mg[0].f[22152]);

    ch->reltol = reltol;

    /* get arrays used in pcg */
    uint8_t* G = mg[0].G;
    double* h = mg[0].h;

    /* reusing mg variables for pcg - q = M\r */
    double* r = mg[0].f;
    double* q = mg[0].v;

    const size_t num_szv = szv[0] * szv[1] * szv[2];
    double* p = (double*)calloc(num_szv, sizeof(double));
    for (int i = 0; i < num_szv; i++)
        p[i] = q[i];


    //mx_zero(q);
    //mx_zero(p);
    for (int i = 0; i < num_szv; i++)
        p[i] = q[i] = 0.0;

    ///* initial guess */
    
    lapmgd(q, x, G, h, mg[0].szv);

    mx_scal(-ONE, q, szv);       /* q = Ax; */
    mx_axpy(-ONE, q, r, szf);    /* r = r - q; */
    res = mx_nrm2(r, szf);

    // printf("INSIDE mx_pcg JUST BEFORE FOR LOOP: res after calling mx_nrm2 = %f\n", res);

    //if (verbose) {
    //    mexPrintf("\niter%6s\tresidual\n", "");
    //}


    // printf("F Matrix inside mgpcg_init  before FOR LOOP...\n");
    // printf("f_matrix at 22149 = %f\n", mg[0].f[22149]);
    // printf("f_matrix at 22150 = %f\n", mg[0].f[22150]);
    // printf("f_matrix at 22151 = %f\n", mg[0].f[22151]);
    // printf("f_matrix at 22152 = %f\n", mg[0].f[22152]);


    printf("=============== STARING FOR LOOP =====================\n");

    for (int32_t i = 0; i < ch->maxit; ++i) {
        if (res <= reltol) {
            ch->isconverged = 1;
            break;
        }

        printf("mx_pcg FOR LOOP at I = %d\n", i);

        prevrho = rho;

        printf("prevrho = %f\n", prevrho);

        for (int i = 0; i < num_szv; i++)
            q[i] = 0.0;

        mx_mg(mg, iter);        /* q = M\r */

        // printf("IN FOR LOOP AFTER mx_mg call : SZ r : [%d] [%d] [%d]\n", mg->szr[0], mg->szr[1], mg->szr[2]);
        // printf("IN FOR LOOP AFTER mx_mg call : SZ q : [%d] [%d] [%d]\n", mg->szr[0], mg->szr[1], mg->szr[2]);

        // Make sure we have size of r and q = num_szv
        rho = mx_dot(r, q, mg[0].szv);     /* rho = r'.q */
        beta = rho / prevrho;

        printf("Updated : rho  = %f\n", rho);
        printf("Updated : beta  = %f\n", beta);

        mx_scal(beta, p, mg[0].szv);
        mx_axpy(ONE, q, p, mg[0].szv);     /* p = q + beta * p */

        lapmgd(q, p, G, h, mg[0].szv);
        mx_scal(-ONE, q, mg[0].szv);       /* q = Ap */
        a = rho / mx_dot(p, q, mg[0].szv); /* a = r'.q / p'.Ap */

        printf("Updated : a  = %f\n", a);

        mx_axpy(a, p, x, mg[0].szv);       /* x = x + a * p */
        mx_axpy(-a, q, r, mg[0].szv);      /* r = r - a * q */

        res = mx_nrm2(r, mg[0].szv);

        printf("Updated : res  = %f\n", res);

        resnorm[i] = res;
        ch->iters += 1;

    //    if (verbose) {
    //        mexPrintf("%5d/%d\t%1.3e\n", i, ch->maxit, res);
    //    }
    }

    //if (verbose) {
    //    mexPrintf("\n");
    //}

    //if (NULL != p) {
    //    mxDestroyArray(p);
    //    p = NULL;
    //}

    r = NULL;
    q = NULL;
    G = NULL;
    h = NULL;

    return;
}

static void
mgpcg_init(struct MG* mg,
    const size_t* szv, const double* mxv,  // V
    const size_t* szf, const double* mxf,  // F
    const size_t* szG, const uint8_t* mxG,  // G
    const size_t* szz, const double* vsz,  // H
    const size_t l)
{
    double* h = NULL;

    size_t size_szv[3] = { szv[0], szv[1], szv[2] };
    
    //mxArray *mxvp = mx_pad_boundary(mxv);
    size_t new_szv[3];
    double* mxvp = mx_pad_boundary(size_szv, mxv, new_szv);

    size_t size_szf[3] = { szf[0], szf[1], szf[2] };
    size_t new_szf[3];
    double* mxfp = mx_pad_boundary(size_szf, mxf, new_szf);

    size_t size_szG[3] = { szG[0], szG[1], szG[2] };
    size_t new_szG[3];
    uint8_t* mxGp = mx_pad_boundary_uint(size_szG, mxG, new_szG);


    //const size_t* szp = (const size_t*)mxGetDimensions(mxvp);
    size_t sz[3] = { new_szv[0], new_szv[1], new_szv[2] };

    double* vsz_dup = (double*)calloc(szz[0]* szz[1]* szz[2], sizeof(double));
    int num_szz = szz[0] * szz[1] * szz[2];
    for (int i = 0; i < num_szz; i++)
        vsz_dup[i] = vsz[i];

    mg[0].v = mxvp;
    mg[0].f = mxfp;
    mg[0].G = mxGp;
    mg[0].h = vsz_dup;
    mg[0].r = (double*)calloc(new_szf[0]* new_szf[1]* new_szf[2], sizeof(double));
    mg[0].B = (uint8_t*)calloc(sz[0]*sz[1]*sz[2], sizeof(uint8_t));

    //Update the sizes after padding
    mg[0].szf[0] = new_szf[0]; mg[0].szf[1] = new_szf[1]; mg[0].szf[2] = new_szf[2];
    mg[0].szv[0] = new_szv[0]; mg[0].szv[1] = new_szv[1]; mg[0].szv[2] = new_szv[2];
    mg[0].szG[0] = new_szG[0]; mg[0].szG[1] = new_szG[1]; mg[0].szG[2] = new_szG[2];
    mg[0].szr[0] = new_szf[0]; mg[0].szr[1] = new_szf[1]; mg[0].szr[2] = new_szf[2];


    boundary_mask(mg[0].B, mg[0].G, sz);


    for (int32_t i = 1; i < l; ++i) {
        sz[0] = ((sz[0] - 2 + 1) >> 1) + 2;
        sz[1] = ((sz[1] - 2 + 1) >> 1) + 2;
        sz[2] = ((sz[2] - 2 + 1) >> 1) + 2;

        mg[i].v = (double*)calloc(sz[0] * sz[1] * sz[2], sizeof(double));
        mg[i].r = (double*)calloc(sz[0] * sz[1] * sz[2], sizeof(double));
        mg[i].f = (double*)calloc(sz[0] * sz[1] * sz[2], sizeof(double));
        mg[i].G = (uint8_t*)calloc(sz[0] * sz[1] * sz[2], sizeof(uint8_t));
        mg[i].B = (uint8_t*)calloc(sz[0] * sz[1] * sz[2], sizeof(uint8_t));

        mg[i].szf[0] = sz[0]; mg[i].szf[1] = sz[1]; mg[i].szf[2] = sz[2];
        mg[i].szv[0] = sz[0]; mg[i].szv[1] = sz[1]; mg[i].szv[2] = sz[2];
        mg[i].szG[0] = sz[0]; mg[i].szG[1] = sz[1]; mg[i].szG[2] = sz[2];
        mg[i].szr[0] = sz[0]; mg[i].szr[1] = sz[1]; mg[i].szr[2] = sz[2];


        coarsen_grid(mg[i].G, mg[i - 1].G, mg[i].szG, mg[i - 1].szG);


        boundary_mask(mg[i].B, mg[i].G, sz);

        mg[i].h = (double*)calloc(num_szz, sizeof(double));
        for (int j = 0; j < num_szz; j++)
            mg[i].h[j] = mg[i - 1].h[j];

        h = mg[i].h;
        h[0] *= 2.0;
        h[1] *= 2.0;
        h[2] *= 2.0;
    }


    h = NULL;

    return;
}

void mgpcg(const uint8_t* szv_in, double* v,
    const uint8_t* szf_in, double* f,
    const uint8_t* szm_in, uint8_t* mask,
    const uint8_t* szz_in, double *vsz,
    double tol,
    int32_t maxit,
    double mgopts_tol,
    int32_t mgopts_maxit,
    int32_t mgopts_mu,
    int32_t mgopts_npre,
    int32_t mgopts_npost,
    int32_t mgopts_nboundary,
    int32_t mgopts_nlevels,
    int32_t in_verbose,
    const uint8_t* szz_out,
    double* out_mat)
{

    //printf("OUT MAT SIZE = [%d %d %d]\n", szz_out[0], szz_out[1], szz_out[2]);

    int MAX_SAMPLES = 20;
    int lin_index = 0;
    int sample_cnt = 0;
    int done = 0;
    for (int i = 0; i < szv_in[0] && !done; i++) {
        for (int j = 0; j < szv_in[1] && !done; j++) {
            for (int k = 0; k < szv_in[2] && !done; k++) {
                lin_index = i * (szv_in[1] * szv_in[2]) + j * szv_in[2] + k;
                if (v[lin_index] > 0) {
                    sample_cnt++;
                    printf("V is non zero at %d\n", lin_index);
                    if (sample_cnt == MAX_SAMPLES) {
                        done = 1;
                    }
                }
            }
        }
    }


    sample_cnt = 0;
    done = 0;
    for (int i = 0; i < szf_in[0] && !done; i++) {
        for (int j = 0; j < szf_in[1] && !done; j++) {
            for (int k = 0; k < szf_in[2] && !done; k++) {
                lin_index = i * (szf_in[1] * szf_in[2]) + j * szf_in[2] + k;
                if (f[lin_index] > 0) {
                    sample_cnt++;
                    //printf("F is non zero at hello %d\n", lin_index);
                    // if (sample_cnt == MAX_SAMPLES) {
                    //     done = 1;
                    // }
                }
            }
        }
    }

    done = 0;
    sample_cnt = 0;
    for (int i = 0; i < szm_in[0] && !done; i++) {
        for (int j = 0; j < szm_in[1] && !done; j++) {
            for (int k = 0; k < szm_in[2] && !done; k++) {
                lin_index = i * (szm_in[1] * szm_in[2]) + j * szm_in[2] + k;
                if (mask[lin_index] > 0) {
                    sample_cnt++;
                    // printf("Mask is non zero at %d\n", lin_index);
                    // if (sample_cnt == MAX_SAMPLES) {
                    //     done = 1;
                    // }
                    
                    // if (sample_cnt < 3*MAX_SAMPLES) {
                    //     printf("Mask is non zero at %d\n", lin_index);
                    // }
                }
            }
        }
    }

    // printf("F Matrix as we get data in...\n");
    // printf("f_matrix at 22149 = %f\n", f[22149]);
    // printf("f_matrix at 22150 = %f\n", f[22150]);
    // printf("f_matrix at 22151 = %f\n", f[22151]);
    // printf("f_matrix at 22152 = %f\n", f[22152]);




    struct PCG ch = {
        .isconverged = 0,
        .iters = 0,
        .resnorm = 0,
        .tol = tol,
        .maxit = maxit
    };

    const size_t szv[] = { szv_in[0], szv_in[1], szv_in[2] };
    const size_t szf[] = { szf_in[0], szf_in[1], szf_in[2] };
    const size_t szm[] = { szm_in[0], szm_in[1], szm_in[2] };
    const size_t szz[] = { szz_in[0], szz_in[1], szz_in[2] };
    
    ch.resnorm = (double*)calloc(ch.maxit, sizeof(double));
    

    const struct Iter iter = {
        .mg = mgopts_maxit,
        .mu = mgopts_mu,
        .pre = mgopts_npre,
        .post = mgopts_npost,
        .bs = mgopts_nboundary,
        .l = mgopts_nlevels,
        .tol = mgopts_tol
    };

    const uint8_t verbose = in_verbose;


    struct MG mg[MAX_DEPTH];
    mgpcg_init(&mg, szv, v, szf, f, szm, mask, szz, vsz, iter.l);



    // printf("F Matrix after mgpcg_init...\n");
    // printf("f_matrix at 22149 = %f\n", mg[0].f[22149]);
    // printf("f_matrix at 22150 = %f\n", mg[0].f[22150]);
    // printf("f_matrix at 22151 = %f\n", mg[0].f[22151]);
    // printf("f_matrix at 22152 = %f\n", mg[0].f[22152]);



    int num_szv = mg[0].szv[0] * mg[0].szv[1] * mg[0].szv[2];
    double *x = (double*)calloc(num_szv, sizeof(double));
    memset(x, 0, num_szv);
    

    //mxArray* x = mxDuplicateArray(mg[0].v); /* initial guess */
    for (int i = 0; i < num_szv; i++)
        x[i] = mg[0].v[i];

    printf("before mpcg\n");

    mx_pcg(x, &ch, mg, &iter, verbose, mg[0].szv, mg[0].szf, mg[0].szG, szz);


    printf("************************ SOLUTION CONVERGED ************************************\n");

    mx_unpad_boundary(x, mg[0].szv, out_mat);

    // Print first 50 non-zero values
    // double *out_v = out_mat;

    // int num_outv = mg[0].szv[0]*mg[0].szv[1]*mg[0].szv[2];
    // printf("NUMBER of output data = %d\n", num_outv);
    // int32_t out_cnt = 0;
    // for (int32_t ii = 0; ii < num_outv; ii++) {
    //     if (out_cnt == 50)
    //         break;
    //     if (out_v[ii] > 0) {
    //         printf("Data = %f\n", out_v[ii]);
    //         out_cnt++;
    //     }
    // }

    //mx_cleanup(x, mg, iter.l);



    return;
}




