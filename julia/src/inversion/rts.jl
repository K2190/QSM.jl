using Images
using MATLAB

include("../utils/fft3.jl")
include("../utils/ifft3.jl")
include("../utils/kernels/dipoleKernel.jl")
include("../utils/kernels/laplaceKernel.jl")
include("../utils/vec.jl")
include("../../third_party/lsmr/lsmr.jl")

mutable struct rts_options
    opts_tolls::Int16
    opts_mu::Float32
    opts_delta::Float32
    opts_tv::Int16
    opts_rtol::Float32
    opts_atol::Float32
    opts_maxit::Int16
    opts_verbose::Int16
    opts_rho::Int16
    opts_alpha::Float32
    opts_rtau::Float32
    opts_stau::Float32
    opts_gamma::Int16
 end     


function rts(f, mask, vsz, bdir)
    #RTS Rapid two-step dipole inversion with sparsity priors.
    #
    #   [x] = RTS(f, mask, vsz, [bdir], [opts]);
    #
    #   Inputs
    #   ------
    #       f               unwrapped local field/phase (3d/4d array).
    #       mask            binary mask of region of interest (3d array).
    #       vsz             voxel size.
    #
    #       bdir            unit vector of B field direction.
    #                       default: [0, 0, 1]
    #       opts.tolls      stopping tolerance, ie # of iterations, for lsmr solver
    #                       default: 4
    #       opts.mu         regularization parameter for tv minimization.
    #                       default: 1e5
    #       opts.delta      threshold for ill-conditioned k-space region.
    #                       default: 0.15
    #       opts.tv         anisotropic (1) or isotropic (2) total variation
    #                       default: 1
    #       opts.rtol       relative stopping tolerance
    #                       default: 1e-1
    #       opts.atol       absolute stopping tolerance
    #                       default: 1e-2
    #       opts.maxit      maximum number of iterations.
    #                       default: 20
    #       opts.verbose    boolean flag for printing information each iteration.
    #                       default: true
    #
    #   Outputs
    #   -------
    #       x               susceptibility map.
    #
    #   Notes
    #   -----
    #       Parameter update and stopping criterion changed from original
    #       publication [1]. Stopping tolerance is now based on primal and dual
    #       residuals instead of the difference of iterates. Parameter update now
    #       includes a clause to keep the primal and dual residuals within a factor
    #       `stau` of each other (in addition to the update scheme from [1]).
    #
    #   References
    #   ----------
    #       [1] Kames C, Wiggermann V, Rauscher A. Rapid two-step dipole inversion
    #       for susceptibility mapping with sparsity priors. Neuroimage. 2018 Feb
    #       15;167:276-83.
    
        # narginchk(3, 5);
    
        # if nargin < 5, opts = []; end
        # if nargin < 4 || isempty(bdir), bdir = [0, 0, 1]; end
    
        # if ~isfield(opts, 'tolls') || isempty(opts.tolls)
        #     opts.tolls = 4;
        # end
        opt_tolls = 4
    
        # if ~isfield(opts, 'mu') || isempty(opts.mu)
        #     opts.mu = 1e5;
        # end
        opt_mu = 1e5
    
        # if ~isfield(opts, 'delta') || isempty(opts.delta)
        #     opts.delta = 0.15;
        # end
        opt_delta = 0.15
    
        # if ~isfield(opts, 'tv') || isempty(opts.tv)
        #     opts.tv = 1; # 1 anisotropic, 2 isotropic
        # end
        opt_tv = 1
    
        # if ~isfield(opts, 'rtol') || isempty(opts.rtol)
        #     opts.rtol = 1e-1;
        # end
        opt_rtol = 1e-1
    
        # if ~isfield(opts, 'atol') || isempty(opts.atol)
        #     opts.atol = 1e-2;
        # end
        opt_atol = 1e-2
    
        # if ~isfield(opts, 'maxit') || isempty(opts.maxit)
        #     opts.maxit = 20;
        # end
        opt_maxit = 20
    
        # if ~isfield(opts, 'verbose') || isempty(opts.verbose)
        #     opts.verbose = 1;
        # end
        opt_verbose = 1
    
        # if ~isfield(opts, 'rho') || isempty(opts.rho)
        #     opts.rho = 10;
        # end
        opt_rho = 10
    
        # # over-relaxation parameter
        # if ~isfield(opts, 'alpha') || isempty(opts.alpha)
        #     opts.alpha = 1.5;
        # end
        opt_alpha = 1.5
    
        # admm parameter update parameters
        # if ~isfield(opts, 'rtau') || isempty(opts.rtau)
        #     opts.rtau = 0.7;
        # end
        opt_rtau = 0.7
    
        # if ~isfield(opts, 'stau') || isempty(opts.stau)
        #     opts.stau = 1e2;
        # end
        opt_stau = 1e2
    
        # if ~isfield(opts, 'gamma') || isempty(opts.gamma)
        #     opts.gamma = 3;
        # end
        opt_gamma = 3

        println("Executing Rapid two-step dipole inversion with sparsity priors...\n")
   
        opts_obj = rts_options(opt_tolls, opt_mu, opt_delta, opt_tv, opt_rtol, opt_atol, 
                              opt_maxit, opt_verbose, opt_rho, opt_alpha, opt_rtau, opt_stau, opt_gamma)
        
    
        # output
        x = zeros_like(size(f), "like", f)

        if measure_time == 1
            t1 = time()
        end
    
        D = dipoleKernel(size(mask), vsz, bdir, "k", eltype(f))

        if measure_time == 1
            t2 = time()
            @printf("====================  Time taken for dipoleKernel = %f msecs\n", (t2-t1)*1000)
        end

        # 4d multi-echo loop
        for t = 1:size(f, 4)
            # if opts.verbose && size(f, 4) > 1
            #     println("Echo: %d/%d\n", t, size(f, 4))
            # end
    
            x[:,:,:,t] = rts_(f[:,:,:,t], mask, D, vsz, opts_obj)
        end

        if measure_time == 1
            tend = time()
            @printf("====================  Total time inside rts == %f msecs\n", (tend-t1)*1000)
        end        
    
       return x
    end
    
    
    function rts_(f, mask, D, vsz, opts)
    
        T = eltype(f)  #class(f)
        sz = size(f)
    
        # lsmr stopping tolerance (iterations)
        tolls = opts.opts_tolls
    
        # stopping tolerance tv
        rtol = opts.opts_rtol
        atol = opts.opts_atol
    
        atolp = sqrt(3*length(f)) * atol
        atold = sqrt(length(f)) * atol
    
        # regularization parameter
        mu = opts.opts_mu
    
        # ill-conditioned k-space threshold
        delta = opts.opts_delta
    
        # admm parameters
        rho = opts.opts_rho
        alpha = opts.opts_alpha
    
        rtau = opts.opts_rtau
        stau = opts.opts_stau
        gamma = opts.opts_gamma
    
        mtv = opts.opts_tv
    
        # admm maxit
        maxit = opts.opts_maxit
    
        # print stuff
        verbose = opts.opts_verbose
    
        # ******************************************************************* #
        # Step 1: Well-conditioned
        # ******************************************************************* #
    
        # reshape to column vector for lsmr
        D = vec_jl(D)

        # tolls iterations of lsmr in kspace
        # x1 = lsmr(@(x,~) D.*x, vec(fft3(f)), 0, eps(T), eps(T), [], tolls, 0, false)
        if measure_time == 1
            t3 = time()
        end


        #Using teh function definition (D.*x) inside lsmr 
        x1 = lsmr(collect(D), collect(vec_jl(fft3(f))), 0, eps(T), eps(T), 1.0e+12, tolls, 0, false)


        if measure_time == 1
            t4 = time()
            @printf("====================  Time taken for LSMR = %f msecs\n", (t4-t3)*1000)
        end

        println("LSMR completed")
    
        x1 = reshape(x1, (sz[1], sz[2], sz[3]))

        x1 = mask .* real(ifft3(x1))

        # reshape back to original size
        D = reshape(D, (sz[1], sz[2], sz[3]))
    
        # ******************************************************************* #
        # Step 2: Ill-conditioned
        # ******************************************************************* #
    
        # mask of well-conditioned frequencies
        M = abs.(D) .> delta
    
        # pre-compute
        irho = 1 / rho
    
        M = mu * M;
        F = M .* fft3(x1)
    
        L = -laplaceKernel(vsz)

        # Try with Padarray and othe image filtering libs, if not lets use MATLAB call
        # L = padarray(L, sz-[3,3,3], "post")
        L = PaddedView(0.0, L, (sz[1], sz[2], sz[3]))
        L = real(fft3(circshift(L, -[1,1,1])));

        iA = 1 ./ (M + rho.*L)

        iA_check = isfinite.(iA)
        iA_check = (!).(iA_check)
        iA[iA_check] .= 0

   
        # initialize admm variables
        x = x1
    
        # p1, p2, p3 = deal(zeros(sz, T))
        
        p1 = zeros(Float64, sz[1], sz[2], sz[3])
        p2 = zeros(Float64, sz[1], sz[2], sz[3])
        p3 = zeros(Float64, sz[1], sz[2], sz[3])

        if T != Float64
            p1 = convert(Array{T}, p1)
            p2 = convert(Array{T}, p2)
            p3 = convert(Array{T}, p3)
        end

        y1, y2, y3 = grad_(x, vsz)
    
        nr = sum(sqrt.(vec(y1.*y1 + y2.*y2 + y3.*y3)))
    
        flag = 0
    
        # if verbose
        #     println("\niter\t\t||r||\t\teps_pri\t\t||s||\t\teps_dual\n");
        # end

        if measure_time == 1
            t5 = time()
            @printf("====================  Time taken prepparing model for iterations = %f msecs\n", (t5-t4)*1000)
        end
    
        for ii = 1:maxit
    
            # *************************************************************** #
            # x - subproblem
            # *************************************************************** #
    
            r1 = rho.*y1 - p1
            r2 = rho.*y2 - p2
            r3 =  rho.*y3 - p3

            r1 = gradAdj_(r1, r2, r3, vsz)
            b  = F + fft3(r1)
            x  = real(ifft3(iA .* b))

            # *************************************************************** #
            # y - subproblem
            # *************************************************************** #
    
            y1_ = y1
            y2_ = y2
            y3_ = y3
    
            dx, dy, dz = grad_(x, vsz)

            # over-relaxation
            dx1 = alpha*dx + (1-alpha)*y1
            dy1 = alpha*dy + (1-alpha)*y2
            dz1 = alpha*dz + (1-alpha)*y3
    
            r1 = dx1 + irho*p1
            r2 = dy1 + irho*p2
            r3 = dz1 + irho*p3

            if mtv == 1
                # anisotropic tv
                y1 = max.(abs.(r1) .- irho, 0) .* sign.(r1)
                y2 = max.(abs.(r2) .- irho, 0) .* sign.(r2)
                y3 = max.(abs.(r3) .- irho, 0) .* sign.(r3)

            else
                # isotropic tv
                s = sqrt(r1.*r1 + r2.*r2 + r3.*r3)
                s = (s - irho) .* (s > irho) ./ s

                s_check = isfinite.(s)
                s_check = (!).(s_check)
                s[s_check] .= 0
    
                y1 = r1 .* s
                y2 = r2 .* s
                y3 = r3 .* s
            end
    
            # *************************************************************** #
            # Lagrange multiplier update
            # *************************************************************** #
    
            p1 = p1 + rho.*(dx1 - y1)
            p2 = p2 + rho.*(dy1 - y2)
            p3 = p3 + rho.*(dz1 - y3)
    
            # *************************************************************** #
            # convergence check
            # *************************************************************** #
    
            nr_ = nr
    
            # primal residual
            r1 = dx - y1
            r2 = dy - y2
            r3 = dz - y3
    
            # dual residual

            s = rho .* gradAdj_(y1-y1_, y2-y2_, y3-y3_, vsz)

            nr = sum(sqrt.(vec(r1.*r1 + r2.*r2 + r3.*r3)))

            # norm(A, p::Real=2)  default is 2
            ns = norm(s)
    
            na = sum(sqrt.(vec(dx.*dx + dy.*dy + dz.*dz)))
            nb = sum(sqrt.(vec(y1.*y1 + y2.*y2 + y3.*y3)))

            ep = atolp + rtol.*max(na, nb)

            # norm(A, p::Real=2)  default is 2
            ed = atold + rtol*norm(rho.*gradAdj_(p1, p2, p3, vsz))
    
            # if verbose
            #     println("%3d/%d\t\t%.4e\t%.4e\t%.4e\t%.4e\n", 
            #         ii, maxit, nr, ep, ns, ed);
            # end
    
            if nr < ep && ns < ed
                break
            end
    
            # *************************************************************** #
            # parameter update
            # *************************************************************** #
    
            if nr > stau*ns || (nr > rtau*nr_ && flag < 3)
                rho  = gamma .* rho;
                irho = 1 / rho;
    
                iA = 1 ./ (M + rho.*L)
                iA_check = isfinite.(iA)
                iA_check = (!).(iA_check)
                iA[iA_check] .= 0
    
                if flag > 5
                    flag = 0;
                end
    
            elseif nr < ns/stau
                rho  = rho ./ gamma^2;
                irho = 1 / rho;
    
                iA = 1 ./ (M + rho.*L);
                iA_check = isfinite.(iA)
                iA_check = (!).(iA_check)
                iA[iA_check] .= 0
    
                flag = flag + 1;
            end
    
        end
    
        # if verbose
        #     fprintf("\n");
        # end
    
        x = mask .* x;
    
        return x
    end
    
    
    function grad_(u, h)
    
        dx = zeros_like(size(u), "like", u)
        dy = zeros_like(size(u), "like", u)
        dz = zeros_like(size(u), "like", u)

        # gradfp_mex(dx, dy, dz, u, h);
        # gradfpf(float *dx, float *dy, float *dz, const float *u, const double *h, const size_t *sz)
        
        hh = zeros(Float64, 3)
        hh[1] = h[1]
        hh[2] = h[2]
        hh[3] = h[3]

        szu = size(u)
        szz = zeros(Int64, 3)
        szz[1] = szu[1]
        szz[2] = szu[2]
        szz[3] = szu[3]


        ccall((:gradfpd, dllpath), Cvoid, (Ref{Float64}, Ref{Float64}, Ref{Float64},
                                           Ref{Float64}, Ref{Float64}, Ref{Int64},),
                                           dx, dy, dz, u, hh, szz)

    
        #ih = 1 ./ h;
    
        #dx = circshift(u, [-1,0,0]) - u;
        #dx = ih(1) .* dx;
    
        #dy = circshift(u, [0,-1,0]) - u;
        #dy = ih(2) .* dy;
    
        #dz = circshift(u, [0,0,-1]) - u;
        #dz = ih(3) .* dz;

        return dx, dy, dz       
    
    end
    
    
    function gradAdj_(dx, dy, dz, h)
    
        d2u = zeros_like(size(dx), "like", dx);
        
        # gradfp_adj_mex(d2u, dx, dy, dz, h);

        hh = zeros(Float64, 3)
        hh[1] = h[1]
        hh[2] = h[2]
        hh[3] = h[3]

        szu = size(d2u)
        szz = zeros(Int64, 3)
        szz[1] = szu[1]
        szz[2] = szu[2]
        szz[3] = szu[3]

        ccall((:gradfp_adjd, dllpath), Cvoid, (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, 
                                               Ref{Float64}, Ref{Int64},), 
                                               d2u, dx, dy, dz, hh, szz)
    

        #ih = -1 ./ h;
    
        #dx = dx - circshift(dx, [1,0,0]);
        #dx = ih(1) .* dx;
    
        #dy = dy - circshift(dy, [0,1,0]);
        #dy = ih(2) .* dy;
    
        #dz = dz - circshift(dz, [0,0,1]);
        #dz = ih(3) .* dz;
    
        #d2u = dx + dy + dz;

        return d2u
    
    end
    