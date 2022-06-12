
include("boundary_mask_l1.jl")
include("check_array.jl")

function erodeMask_(mask, vx)

    # try
        # fast implementation in C: ./poisson_solver/mex/boundary_mask_l1_mex.c
        #b = false(size(mask));
        b = zeros(size(mask))
        bUint = convert(Array{UInt8}, b)

        mask = convert(Array{UInt8}, mask)
        
        szu = size(b)
        szz = zeros(Int64, 3)
        szz[1] = szu[1]
        szz[2] = szu[2]
        szz[3] = szu[3]

        ccall((:boundary_mask_l1, dllpath), Cvoid, (Ref{UInt8}, Ref{UInt8}, Ref{Int64}, UInt8,),
                                                                          bUint, mask, szz, vx)


        #<TODO> Enable c code
       # boundary_mask_l1_mex(b, mask, vx);

       mask = mask - bUint;
       return mask

    # catch ex
    #     # println(ex.identifier, "WARNING %s\n", ex.message, "using fallback method");
    #     # for ii = 1:vx
    #     #     mask = mask - bwperim(mask, 6);
    #     #     return mask
    #     # end
    #     return mask
    # end

    return 0
end

function erodeMask(mask, vx)
#ERODEMASK Erode 3d binary mask.
#
#   [mask] = ERODEMASK(mask, vx);
#
#   Inputs
#   ------
#       mask    3d binary mask.
#       vx      number of voxels to erode.
#
#   Outputs
#   -------
#       mask    eroded binary mask.
#
#   See also DILATEMASK

    #<TODO> Find Julia equivalent narginchk
    #narginchk(2, 2)

    validateinputs(mask, vx)

    if vx < 1
        return
    end

    #T = class(mask);
    T = typeof(mask)
    #mask = logical(mask);

    mask = mask .!= 0

    mask = erodeMask_(mask, vx);

    #mask = cast(mask, T);
    mask = convert(T, mask)

    return mask

end





function validateinputs(mask, vx)

    # classes = {"logical", "numeric"};
    # attributes = {"real", "ndims", 3, "finite", "binary"};
    # validateattributes(mask, classes, attributes, mfilename, "mask", 1);    

    # classes = {"numeric"};
    # attributes = {"real", "scalar", "finite", "nonnegative", "integer"};
    # validateattributes(vx, classes, attributes, mfilename, "vx", 2);

end
