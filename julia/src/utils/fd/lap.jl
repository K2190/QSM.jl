
function lap_(u, mask, vsz)

    #d2u = zeros(size(u), "like", u);
    d2u = zeros(size(u))    # <TODO> check if "like" is needed
    lap1_mex(d2u, u, mask, vsz);

    return d2c
end

function [] = validateinputs(u, mask, vsz)

    classes = {"single", "double"};
    attributes = {"real", "ndims", 3, "finite"};
    validateattributes(u, classes, attributes, mfilename, 'u', 1);

    classes = {"logical", "numeric"};
    attributes = {"real", "ndims", 3, "size", size(u), "finite", "binary"};
    validateattributes(mask, classes, attributes, mfilename, "mask", 2);

    classes = {"numeric"};
    attributes = {"real", "vector", "numel", 3, "finite", ">", 0};
    validateattributes(vsz, classes, attributes, mfilename, "vsz", 3);

end

function lap(u, mask, vsz)
#LAP Second order central difference Laplacian.

#   [d2u] = LAP(u, [mask], [vsz]);

    #<TODO> Find Julia equivalent narginchk
    #narginchk(1, 3);

    #<TODO> Pass these as default argiments from caller
    #if nargin < 3,  vsz = [1, 1, 1]; end
    #if nargin < 2 || isempty(mask), mask = true(size(u(:,:,:,1))); end

    validateinputs(u, mask, vsz)

    #mask = logical(mask)
    mask = !=(0).(mask)

    #vsz = double(vsz)
    vsz = convert(Float64, vsz)


    d2u = lap_(u, mask, vsz)

    return d2u

end






