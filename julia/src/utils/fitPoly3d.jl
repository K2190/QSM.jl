using MATLAB

include("vec.jl")
include("nd_grid.jl")
include("write_variable.jl")

function find_E_qr(inv_A, Q, R)
    P1 = Q * R;
    P2 = inv_A * P1
    P2 = round.(P2)

    num = rank(Q)
    E = zeros(Int32, num)

    for i = 1:num
        P2_2 = P2[:, i]
        indx = findall(x->x==1, P2_2)
        E[i] = indx
    end
end

function calculateQER(A)
    

    szzA = size(A)
    newSzA = 0
    if szzA[1] > szzA[2]
        newSzA = szzA[1]
    else
        newSzA = szzA[2]
    end

    Anew = zeros(Float64, newSzA, newSzA)
    #Fist copy original matrix
    for i = 1:szzA[1]
        for j=1:szzA[2]
            Anew[i, j] = A[i, j]
        end
    end
    
    if szzA[1] > szzA[2]
        for i = 1:szzA[1]  
            for j = szzA[2]+1:szzA[1]
                Anew[i, j] = 0.0
            end
        end
    else
        for i = 1:szzA[2]  
            for j = szzA[1]+1:szzA[1]
                Anew[i, j] = 0.0
            end
        end
    end

    invA = inv(Anew)

    Q, R = qr!(Anew, Val(true))
    E = find_E_qr(invA, Q, R)

    return [Q, R, E]
end


function mergeArraysHorizontally(arr1, arr2)
    sz1 = size(arr1)
    sz2 = size(arr2)
    numr = sz1[1]
    numc = sz1[2] + sz2[2]
    m = zeros(Float64, numr, numc)
    for i = 1 : numr
        for j = 1 : sz1[2]
            m[i, j] = arr1[i, j]
        end
    end
    for i = 1 : numr
        for j = 1:sz2[2]
            m[i, j+sz1[2]] = arr2[i, j]
        end
    end
    return m
end


function buildmodel(order, p)
    if p == 0
        m = []
  
    elseif order == 0
        m = zeros(Float64, 1, p)
  
    elseif p == 1
      szm = order + 1
      m = zeros(Float64, szm, 1)
      for i = 1 : szm
          m[i, 1] = szm - i
      end
  
    else
        m = Array{Float64}(undef, 0, p)
        for k = order:-1:0
            t = buildmodel(order-k, p-1)
            X = [k]
            sz1 = size(t, 1)
            temp = repeat(X, outer = [sz1, 1])
            temp2 = mergeArraysHorizontally(temp, t)
            m = [m; temp2]
        end
    end
  
    return m
  end

function validateinputs(x, n, mask, vsz)

    # sz = size(x);

    # if ndims(x) < 4
    #     nd = 3;
    # else
    #     nd = 4;
    # end

    # classes = {"single", "double"};
    # attributes = {"real", "ndims", nd, "finite"};
    # validateattributes(x, classes, attributes, mfilename, "x", 1);

    # classes = {"numeric"};
    # attributes = {"real", "scalar", "finite", "integer", "nonnegative"};
    # validateattributes(n, classes, attributes, mfilename, "n", 2);

    # classes = {"logical", "numeric"};
    # attributes = {"real", "ndims", 3, "size", sz[1:3], "finite", "binary"};
    # validateattributes(mask, classes, attributes, mfilename, "mask", 3);

    # classes = {"numeric"};
    # attributes = {"real", "vector", "numel", 3, "finite", ">", 0};
    # validateattributes(vsz, classes, attributes, mfilename, "vsz", 4);

end

function get_diagonal(inp_mat)
    szz =size(inp_mat)
    num = szz[1]
    diag_mat = zeros(eltype(inp_mat), num)
 
    for i=1:num
     diag_mat[i] = inp_mat[i,i]
    end
 
    return diag_mat
 end


function fitPoly3d(x, n, mask, vsz)
#FITPOLY3D Fit 3d polynomial of order n to 3d multi-echo data
#
#   [y] = FITPOLY3D(x, n, [mask], [vsz]);
#
#   Inputs
#   ------
#       x       multi-echo data, 3d/4d array.
#       n       order of polynomial to fit.
#
#       mask    binary mask of voxels to include in fit.
#               default = ones(size(x))
#       vsz     voxel size.
#               default = [1, 1, 1]
#
#   Outputs
#   -------
#       y       polynomial fit of x
#
#   Notes
#   -----
#       Adapted from John D'Errico's polyfitn toolbox:
#       https://www.mathworks.com/matlabcentral/fileexchange/34765-polyfitn
#
#   See also POLYFITN, FITPOLY2D

    #<TODO> Find Julia equivalent narginchk
    #narginchk(2, 4);

    #<TODO> Pass these as default argiments from caller
    #if nargin < 4, vsz = [1, 1, 1]; end
    #if nargin < 3 || isempty(mask), mask = true(size(x(:,:,:,1))); end

    validateinputs(x, n, mask, vsz)

    if measure_time == 1
        t0 = time()
    end


    #T = class(x);
    T = typeof(x[1])

    sz = size(mask);

    #sz = cast(sz, T);   MATLAB
    # sz = convert(T, sz)  JULIA, may not be needed now

    #vsz = cast(vsz, T);
    vsz = convert(Array{T}, vsz)

    #mask = logical(mask);

    mask = mask .!= 0

    # build the model
    k = buildmodel(n, 3);
    #k = cast(k, T);  MATLAB
    #k = convert(T, k)  JULIA, may not be needed now

    # generate grid
    sz = collect(sz)
    low = floor.(sz / 2);
    #high = low - ~mod(sz, 2);
    temp = mod.(sz,2)
    temp = temp == 1 ? 0 : 1;
    high = low .- temp

    low = vsz .* low;
    high = vsz .* high;

    xx = collect(-low[1] : vsz[1] : high[1] )
    yy = collect(-low[2] : vsz[2] : high[2] )
    zz = collect(-low[3] : vsz[3] : high[3] )

    if measure_time == 1
        t1 = time()
        @printf("====================  Time taken for buildmodel = %f msecs\n", (t1-t0)*1000)
     end

    # xx = range(-low[1] , high[1], step=vsz[1]) |> collect
    # yy = range(-low[2] , high[2], step=vsz[2]) |> collect
    # zz = range(-low[3] , high[3], step=vsz[3]) |> collect
    X, Y, Z = nd_grid(xx, yy, zz)

    if measure_time == 1
        t2 = time()
        @printf("====================  Time taken for nd_grid = %f msecs\n", (t2-t1)*1000)
    end

    # X, Y, Z = ndgrid(-low[1]:vsz[1]:high[1], 
    #                    -low[2]:vsz[2]:high[2], 
    #                    -low[3]:vsz[3]:high[3])


    # scale to unit variance
    I = [X[mask] Y[mask] Z[mask]];


    s = sqrt.(get_diagonal(cov(I)))

    C = findall(s .== 0)
    s[C] .= 1
    diagS = Diagonal(1 ./ s)
    I = I * diagS


    if measure_time == 1
        t3 = time()
        @printf("====================  Time taken for Diagonal matrix = %f msecs\n", (t3-t3)*1000)
    end

    # build the design matrix
    A = ones(T, size(I, 1), size(k, 1));

    for ii = 1:size(A, 2)
        A[:,ii] = I[:,1].^k[ii, 1] .* I[:,2].^k[ii, 2] .* I[:,3].^k[ii, 3]
    end


    if measure_time == 1
        t4 = time()
        @printf("====================  Time taken for design matrix = %f msecs\n", (t4-t3)*1000)
    end


    println("Generating QR...\n")

    # Using MATLAB direct call
    Q, R, E = mxcall(:qr, 3, A, 0)
    E = convert(Array{Int8}, E)


    if measure_time == 1
        t5 = time()
        @printf("====================  Time taken for QR = %f msecs\n", (t5-t4)*1000)
    end

    I = [vec_jl(X) vec_jl(Y) vec_jl(Z)];
    diagS2 = Diagonal(1 ./ s)
    I = I * diagS2

    # build the design matrix
    A = ones(T, size(I, 1), size(k, 1));

    for ii = 1:size(A, 2)
        A[:,ii] = I[:,1].^k[ii, 1] .* I[:,2].^k[ii, 2] .* I[:,3].^k[ii, 3]
    end


    if measure_time == 1
        t6 = time()
        @printf("====================  Time taken for design matrix = %f msecs\n", (t6-t5)*1000)
    end
 
    y = zeros(T, sz[1], sz[2], sz[3]);

    println("Generating poly fit...\n")

    for t = size(x, 4):-1:1
        x_ = x[:,:,:,t]
        xm = x_[mask]
        xm = vec_jl(xm)
        Qxm = transpose(Q)*xm
        temp22 = R \ Qxm
        c = rand(1,35)
        c[E] = temp22
        y_ = A * vec_jl(c)
        y[:,:,:,t] = reshape(y_, sz[1], sz[2], sz[3])
    end


    if measure_time == 1
        t7 = time()
        @printf("====================  Time taken for Generating poly fit = %f msecs\n", (t7-t6)*1000)
    end


    if measure_time == 1
        @printf("====================  Total time inside poly fit = %f msecs\n", (t7-t0)*1000)
    end

    return y
end

