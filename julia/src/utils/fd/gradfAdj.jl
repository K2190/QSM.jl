
function gradAdj_(x, y, z, mask, vsz)

    #du = zeros(size(x), "like", x);
    du = zeros(size(x));

    if all(vec(mask))
        gradf_adj_mex(du, x, y, z, vsz);
    else
        gradfm_adj_mex(du, x, y, z, mask, vsz);
    end

    return du
end


function gradAdj_(x, y, z, h)

    ih = -1 ./ h

    d = x - Base.circshift(x, [1,0,0]);
    d[1,:,:] = d[2,:,:]
    d = ih[1] .* d;

    du = d;

    d = y - Base.circshift(y, [0,1,0]);
    d[:,1,:] = d[:,2,:]
    d = ih[2] .* d;

    du = du + d;

    d = z - Base.circshift(z, [0,0,1]);
    d[:,:,1] = d[:,:,2]
    d = ih[3] .* d;

    du = du + d;

    return du
end



function [] = validateinputs(x, y, z, mask, vsz)

    sz = size(x);

    classes = {"single", "double"};
    attributes = {"real", "ndims", 3, "finite"};
    validateattributes(x, classes, attributes, mfilename, "x", 1);

    classes = {"single", "double"};
    attributes = {"real", "ndims", 3, "size", sz, "finite"};
    validateattributes(y, classes, attributes, mfilename, "y", 2);

    classes = {"single", "double"};
    attributes = {"real", "ndims", 3, "size", sz, "finite"};
    validateattributes(z, classes, attributes, mfilename, "z", 3);

    classes = {"logical", "numeric"};
    attributes = {"real", "ndims", 3, "size", sz, "finite", "binary"};
    validateattributes(mask, classes, attributes, mfilename, "mask", 4);

    classes = {"numeric"};
    attributes = {"real", "vector", "numel", 3, "finite", ">", 0};
    validateattributes(vsz, classes, attributes, mfilename, "vsz", 5);

end



function [du] = gradfAdj(x, y, z, mask, vsz)
#GRADFADJ Adjoint of first order forward difference gradient.
#
#   [du] = GRADFADJ(x, y, z, [mask], [vsz]);
#
#   See also GRADF

    #<TODO> Find Julia equivalent narginchk
    #narginchk(3, 5);

    #<TODO> Pass these as default argiments from caller
    #if nargin < 5,  vsz = [1, 1, 1]; end
    #if nargin < 4 || isempty(mask), mask = true(size(x(:,:,:,1))); end

    validateinputs(x, y, z, mask, vsz);

    # Single-precision variables in MATLAB® are stored as 4-byte (32-bit) 
    # floating-point values
    #if isa(x, "single") || isa(y, "single") || isa(z, "single")
    if check_precision(x, "single") || check_precision(y, "single") || check_precision(z, "single")
        #x = cast(x, "single");
        #y = cast(y, "single");
        #z = cast(z, "single");
        if isinteger(x)
            x = convert(Int32, x)
        else
            x = convert(Float32, x)
        if isinteger(y)
            y = convert(Int32, y)
        else
            y = convert(Float32, y)
        if isinteger(z)
            z = convert(Int32, z)
        else
            z = convert(Float32, z)
    end

    #mask = logical(mask);
    mask = !=(0).(mask)

    #vsz = double(vsz);
    vsz = convert(Float64, vsz)

    du = gradAdj_(x, y, z, mask, vsz);

end

