#include <inttypes.h>
#include <string.h>
#include <stddef.h>
#include "boundary_mask_mex.h"
#include "coarsen_grid_mex.h"
#include "correct_mex.h"
#include "gauss_seidel_mex.h"
#include "mx_blas.h"
#include "mx_util.h"
#include "prolong_mex.h"
#include "restrict_mex.h"
#include "residual_mex.h"


#define MAX_DEPTH 16


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

