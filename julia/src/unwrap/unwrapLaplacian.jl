
include("../utils/matlab_2_julia.jl")
include("../utils/cropIndices.jl")
include("lapw.jl")
include("../utils/poisson_solver/mgpcg.jl")


mutable struct mgoptions
    maxit::UInt16
    npre::UInt16
    npost::UInt16
    nboundary::UInt16
end

function unwrapLaplacian(phas, mask, vsz)
#UNWRAPLAPLACIAN Laplacian phase unwrapping.
#
#   [uphas, d2phas] = UNWRAPLAPLACIAN(phas, mask, vsz);
#
#   Inputs
#   ------
#       phas    wrapped phase (3d/4d array).
#       mask    binary mask of region of interest (3d array).
#       vsz     voxel size.
#
#   Outputs
#   -------
#       uphas   unwrapped local phase.
#       d2phas  Laplacian of unwrapped phase.
#
#   Notes
#   -----
#       The Laplacian is computed using second order central finite differences
#       on the complex phase (see lapw.m).
#
#       The resulting Poisson's equation is then solved inside the ROI with
#       homogenous Dirichlet BCs (uphas(~mask) = 0) [2] using a multigrid
#       preconditioned conjugate gradient method (mgpcg.m). The boundary of the
#       ROI is set such that values outside of it (mask = 0) are taken as
#       boundary points and values inside of it (mask = 1) as interior points.
#
#       This method combines phase unwrapping [1] and harmonic background field
#       removing [2] (see lbv.m).
#
#   References
#   ----------
#       [1] Schofield MA, Zhu Y. Fast phase unwrapping algorithm for
#       interferometric applications. Optics letters. 2003 Jul 15;28(14):1194-6.
#
#       [2] Zhou D, Liu T, Spincemaille P, Wang Y. Background field removal by
#       solving the Laplacian boundary value problem. NMR in Biomedicine. 2014
#       Mar;27(3):312-9.
#
#   See also LAPW, MGPCG

    #narginchk(3, 3)

    # validateinputs(phas, mask, vsz)

    #mask = logical(mask)

    #<TODO> REMOVE LATER, disabled for DEBUGGING
    #mask = !=(0).(mask)

    vsz = convert_to_double(vsz)

    # output
    uphas = zeros_like(size(phas), "like", phas)

    # <TODO> Check what is this after debugging
    # if nargout > 1
    #     d2phas = zeros_like(size(phas), "like", phas)
    # end
    d2phas = zeros_like(size(phas), "like", phas)

    # crop to avoid unnecessary work. phas(~mask) = 0
    ix, iy, iz = cropIndices(mask)

    # ix = vec(max(1, ix(1)-1) : min(size(mask, 1), ix(end)+1));

    min_1 = min(size(mask, 1), ix[end]+1)
    max_1 = max(1, ix[1]-1)
    ix = collect(max_1:min_1)

    min_1 = min(size(mask, 2), iy[end]+1)
    max_1 = max(1, iy[1]-1)
    iy = collect(max_1:min_1)

    min_1 = min(size(mask, 3), iz[end]+1)
    max_1 = max(1, iz[1]-1)
    iz = collect(max_1:min_1)

    mask = mask[ix, iy, iz]

    for t = 1:size(phas, 4)
        uphas[ix,iy,iz,t], d2p = unwrap_(phas[ix,iy,iz,t], mask, vsz)

        # if nargout > 1
        #     d2phas[ix,iy,iz,t] = d2p
        # end
        d2phas[ix,iy,iz,t] = d2p
    end

    return uphas, d2phas
end


function unwrap_(u, m, h)

    # get Laplacian of unwrapped phase from wrapped phase
    d2u = lapw(u, h);

    # solve Poisson's equation with homogenous Dirichlet BCs:
    #   -Delta u = f,   for x in mask
    #   u(~mask) = 0

    # TODO: input params
    tolcg = 1e-8
    maxitcg = ceil(sqrt(length(m)))

    mgopts_maxit = 2
    mgopts_npre = 2
    mgopts_npost = 2
    mgopts_nboundary = 2

    mgopts_obj = mgoptions(mgopts_maxit, mgopts_npre, mgopts_npost, mgopts_nboundary)

    # u = mgpcg(-m.*d2u, m, h, tolcg, maxitcg, mgopts);  This is in matlab
    # u = mgpcg(-m.*d2u, m, h, tolcg, maxitcg, mgopts_obj, true)
    # Try this  d2u_val = d2u[1]....m.*d2u_val

    u = mgpcg(-m.*d2u, m, h, tolcg, maxitcg, mgopts_obj, true)

    return u, d2u
end


function validateinputs(phas, mask, vsz)

    # sz = size(phas)

    # if ndims(phas) < 4
    #     nd = 3
    # else
    #     nd = 4
    # end

    # classes = {"single", "double"}
    # attributes = {"real", "ndims", nd, "finite"}
    # validateattributes(phas, classes, attributes, mfilename, "phas", 1)

    # classes = {"logical", "numeric"}
    # attributes = {"real", "ndims", 3, "size", sz(1:3), "finite", "binary"}
    # validateattributes(mask, classes, attributes, mfilename, "mask", 2)

    # classes = {"numeric"}
    # attributes = {"real", "vector", "numel", 3, "finite", ">", 0}
    # validateattributes(vsz, classes, attributes, mfilename, "vsz", 3)

end
