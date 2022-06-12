include("../write_variable.jl")

function mgpcg(f, mask, vsz, tol_in, maxit, mgopts, verbose)
#MGPCG Multigrid preconditioned conjugate gradient for Poisson's equation
#   with homogenous Dirichlet BCs.
#
#                           -Delta u = f, for x in mask
#                           u(~mask) = 0
#
#   [v] = MGPCG(f, mask, vsz, [tol], [maxit], [mgopts], [verbose]);
#
#   Inputs
#   ------
#       f                   right hand side (3d array): -Delta u = f.
#       mask                binary mask (3d array). 1 = interior, 0 = Dirichlet.
#       vsz                 voxel size.
#
#       tol                 stopping tolerance for conjugate gradient method.
#                           ||A*v_k - f||_2 <= tol * ||A*v_0 - f||_2
#                           default = sqrt(eps(class(f)))
#       maxit               maximum number of conjugate gradient iterations.
#                           default = numel(f)
#       mgopts.tol          stopping tolerance for multigrid cycles.
#                           ||A*v_k - f||_2 <= tol * ||A*v_0 - f||_2
#                           default = -1 (fixed iterations)
#       mgopts.maxit        maximum number of multigrid cycles.
#                           default = 1
#       mgopts.mu           type of mu-cycle: 1 = V-cycle,  2 = W-cycle, ...
#                           default = 1
#       mgopts.npre         number of pre-relaxation sweeps.
#                           default = 1
#       mgopts.npost        number of post-relaxation sweeps.
#                           default = 1
#       mgopts.nboundary    number of extra boundary sweeps after interior
#                           sweeps on the downstroke and before interior sweeps
#                           on the upstroke.
#                           default = 2
#       mgopts.nlevels      number of multigrid levels.
#                           default = max(1, floor(log2(min(size(v)))) - 3)
#
#   Notes
#   -----
#       The Poisson equation is discretized using a standard second order
#       central finite-difference stencil.
#
#       The default smoother is a parallel red-black Gauss-Seidel. A serial
#       Gauss-Seidel is also implemented and can be used by setting the
#       gs_red_black flag in make.m to 0 and recompiling.
#
#       Coarse grids are generated such that a coarse grid point is a Dirichlet
#       point if any of its eight fine children is a Dirichlet point. The
#       coarse grid point is an interior point otherwise.
#
#       Prolongation and restriction operators are trilinear interpolation and
#       full-weighting (27 points), respectively. Prolongation and restriction
#       operators only prolong and restrict into interior grid points, a value
#       of 0 is assigned otherwise.
#
#       Inputs are automatically padded with a boundary layer of Dirichlet
#       points.
#
#   References
#   ----------
#       [1] Briggs, W. L., & McCormick, S. F. (2000). A multigrid tutorial
#       (Vol. 72). Siam.
#
#   See also CG, MG, FMG

    # narginchk(3, 7);

    # if nargin < 7, verbose = false; end
    # if nargin < 6 || isempty(mgopts), mgopts = []; end
    # if nargin < 5 || isempty(maxit), maxit = numel(f); end
    # if nargin < 4 || isempty(tol), tol = sqrt(eps(class(f))); end

    # if isempty(mgopts.tol)
    #     mgopts.tol = -1 # fixed iterations
    # end

    # if isempty(mgopts.maxit)
    #     mgopts.maxit = 1
    # end

    # if isempty(mgopts.mu)
    #     mgopts.mu = 1 # 1 = V-cycle, 2 = W-cycle
    # end

    # if isempty(mgopts.npre)
    #     mgopts.npre = 1
    # end

    # if isempty(mgopts.npost)
    #     mgopts.npost = 1
    # end

    # if isempty(mgopts.nboundary)
    #     mgopts.nboundary = 2
    # end

    # if isempty(mgopts.nlevels)
    #     mgopts.nlevels = max(1, floor(log2(min(size(f)))) - 3)
    # end


    # validateinputs(f, mask, vsz, tol, maxit, mgopts, verbose)

    #mask = logical(mask)
    #mask = !=(0).(mask)

    vsz = convert_to_double(vsz)

    # initial guess
    v = zeros_like(size(f), "like", f)

    println("********** Executing MGPCG ************")

    # We need the size inside
    szv_ = size(v)
    szf_ = size(f)
    szm_ = size(mask)

    vsz = [ 1.0, 1.0, 1.0 ]
    szv = [ szv_[1], szv_[2], szv_[3] ]
    szf = [ szf_[1], szf_[2], szf_[3] ]
    szm = [ szm_[1], szm_[2], szm_[3] ]
    szz = [ 3, 1, 1 ]

    tol              = tol_in
    maxit            = maxit
    mgopts_maxit     = mgopts.maxit
    mgopts_npre      = mgopts.npre
    mgopts_npost     = mgopts.npost
    mgopts_mu        = 1
    mgopts_nboundary = mgopts.nboundary
    mgopts_nlevels   = 2
    mgopts_tol       = -1.0
    verbose          = 0


    szv              = convert(Array{UInt8}, szv)
    szf              = convert(Array{UInt8}, szf)
    szm              = convert(Array{UInt8}, szm)
    szz              = convert(Array{UInt8}, szz)
    tol              = convert(Float64, tol)
    mgopts_tol       = convert(Float64, mgopts_tol)
    max_it           = convert(Int32, maxit)
    mgopts_maxit     = convert(Int32, mgopts_maxit)
    mgopts_npre      = convert(Int32, mgopts_npre)
    mgopts_npost     = convert(Int32, mgopts_npost)
    mgopts_mu        = convert(Int32, mgopts_mu)
    mgopts_nboundary = convert(Int32, mgopts_nboundary)
    mgopts_nlevels   = convert(Int32, mgopts_nlevels)
    verbose          = convert(Int32, verbose)


    sz_outv  = [ szv_[1], szv_[2], szv_[3] ]
    sz_outv  = convert(Vector{UInt8}, sz_outv)
    out_v = zeros(Float64, szv_[1]*szv_[2]*szv_[3])

    mask = convert(Array{UInt8}, mask)
    v = convert(Array{Float64}, v)
    f = convert(Array{Float64}, f)

    ccall((:mgpcg, dllpath), Ref{Float64},  (Ref{UInt8}, Ref{Float64}, 
                                             Ref{UInt8}, Ref{Float64},
                                             Ref{UInt8}, Ref{UInt8},
                                             Ref{UInt8}, Ref{Float64},
                                             Float64, Int32,
                                             Float64,
                                             Int32,
                                             Int32,
                                             Int32,
                                             Int32,
                                             Int32,
                                             Int32,
                                             Int32,
                                             Ref{UInt8}, Ref{Float64},), 
                                             szv, v, 
                                             szf, f,
                                             szm, mask,
                                             szz, vsz,
                                             tol,
                                             max_it,
                                             mgopts_tol,
                                             mgopts_maxit,
                                             mgopts_mu,
                                             mgopts_npre,
                                             mgopts_npost,
                                             mgopts_nboundary,
                                             mgopts_nlevels,
                                             verbose,
                                             sz_outv, out_v)


    out_v = reshape(out_v, (szv_[1], szv_[2], szv_[3]))

    return out_v
end 


function validateinputs(f, mask, vsz, tol, maxit, o, verbose)

    # sz = size(f)

    # classes = {"single", "double"}
    # attributes = {"real", "ndims", 3, "size", sz, "finite"}
    # validateattributes(f, classes, attributes, mfilename, "f", 1)

    # classes = {"logical", "numeric"}
    # attributes = {"real", "ndims", 3, "size", sz, "binary"}
    # validateattributes(mask, classes, attributes, mfilename, "mask", 2)

    # classes = {"numeric"}
    # attributes = {"real", "vector", "numel", 3, "finite", ">", 0}
    # validateattributes(vsz, classes, attributes, mfilename, "vsz", 3)

    # # cg opts
    # classes = {"numeric"}
    # attributes = {"real", "scalar", "finite"}
    # validateattributes(tol, classes, attributes, mfilename, "tol", 4)

    # classes = {"numeric"}
    # attributes = {"scalar", "integer", ">", 0}
    # validateattributes(maxit, classes, attributes, mfilename, "maxit", 5)

    # classes = {"logical", "numeric"}
    # attributes = {"real", "scalar", "binary"}
    # validateattributes(verbose, classes, attributes, mfilename, "verbose", 7)

    # mg opts
    # classes = {"numeric"}
    # attributes = {"real", "scalar", "finite"}
    # validateattributes(o.tol, classes, attributes, 
    #     mfilename, "opts.tol", 6)

    # classes = {"numeric"}
    # attributes = {"scalar", "integer", ">", 0}
    # validateattributes(o.maxit, classes, attributes, 
    #     mfilename, "opts.maxit", 6)

    # classes = {"numeric"}
    # attributes = {"scalar", "integer", ">", 0}
    # validateattributes(o.mu, classes, attributes, 
    #     mfilename, "opts.mu", 6)

    # classes = {"numeric"}
    # attributes = {"scalar", "integer", ">=", 0}
    # validateattributes(o.npre, classes, attributes, 
    #     mfilename, "opts.npre", 6)

    # classes = {"numeric"}
    # attributes = {"scalar", "integer", ">=", 0}
    # validateattributes(o.npost, classes, attributes, 
    #     mfilename, "opts.npost", 6)

    # classes = {"numeric"}
    # attributes = {"scalar", "integer", ">=", 0}
    # validateattributes(o.nboundary, classes, attributes, 
    #     mfilename, "opts.nboundary", 6)

    # classes = {"numeric"}
    # attributes = {"scalar", "integer", ">", 0, "<", floor(log2(min(sz)))-1}
    # validateattributes(o.nlevels, classes, attributes, 
    #     mfilename, "opts.nlevels", 6);

end
