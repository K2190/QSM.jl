function lapw(phas, vsz)
#LAPW Discrete Laplacian of wrapped phase data.
#
#   [d2phas] = LAPW(phas, vsz);
#
#   Inputs
#   ------
#       phas    wrapped phase (3d array).
#       vsz     voxel size (3 element vector).
#
#   Outputs
#   -------
#       d2phas  Laplacian of unwrapped phase.
#
#   Notes
#   -----
#       The Laplacian is computed using second order central finite differences
#       on the complex phase:
#
#           d2u/dx2 = arg(exp( i * [u(x-1) - 2u(x) + u(x+1)] ))

    #narginchk(2, 2);

    validateinputs(phas, vsz)
    vsz = convert_to_double(vsz)

    d2phas = lapw_(phas, vsz)

    return d2phas
end


function lapw_(u, h)

    d2u = zeros_like(size(u), "like", u)

    T1 = typeof(u)
    T2 = typeof(h)

    hh = zeros(Float64, 3)
    hh[1] = h[1]
    hh[2] = h[2]
    hh[3] = h[3]

    szu = size(u)
    szz = zeros(Int64, 3)
    szz[1] = szu[1]
    szz[2] = szu[2]
    szz[3] = szu[3]



    #lapw_mex(d2u, u, h);
    ccall((:lapwd, dllpath), Cvoid, (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64},),
                                     d2u, u, hh, szz)
    
    

    return d2u
end




function validateinputs(phas, vsz)

    # classes = {"single", "double"};
    # attributes = {"real", "ndims", 3, "finite"};
    # validateattributes(phas, classes, attributes, mfilename, "phas", 1);

    # classes = {"numeric"};
    # attributes = {"real", "vector", "numel", 3, "finite", ">", 0};
    # validateattributes(vsz, classes, attributes, mfilename, "vsz", 2);

end
