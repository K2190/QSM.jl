function gradf(u, mask, vsz)
#GRADF First order forward difference gradient.
#   Boundary condition: u^{k+1} - 2u^{k} + u^{k-1} = 0.
#
#   [dx, dy, dz] = GRADF(u, [mask], [vsz]);
#
#   See also GRADB, GRADC, GRADFADJ

    #<TODO> Find Julia equivalent narginchk
    #narginchk(1, 3);

    #<TODO> Pass these as default argiments from caller
    #if nargin < 3,  vsz = [1, 1, 1]; end
    #if nargin < 2 || isempty(mask), mask = true(size(u)); end

    validateinputs(u, mask, vsz);

    #mask = logical(mask);
    mask = !=(0).(mask)

    #vsz = double(vsz);
    vsz = convert(Float64, vsz)

    dx, dy, dz = grad_(u, mask, vsz);

    return [dx, dy, dz]
end


function grad_(u, mask, vsz)

    #dx = zeros(size(u), 'like', u);
    #dy = zeros(size(u), 'like', u);
    #dz = zeros(size(u), 'like', u);

    dx = zeros(size(u))
    dy = zeros(size(u))   
    dz = zeros(size(u)) 
    dx = convert(typeof(u), dx)
    dy = convert(typeof(u), dy)
    dz = convert(typeof(u), dz)

    #if all(vec(mask))
    if all(!=(0), mask)
        gradf_mex(dx, dy, dz, u, vsz);
    else
        gradfm_mex(dx, dy, dz, u, mask, vsz);
    end

    return [dx, dy, dz]
end


function grad_(u, h)

    ih = 1 ./ h;

    dx = Base.circshift(u, [-1, 0, 0]) - u
    dx[end,:,:] = dx[end-1,:,:]
    dx = ih[1] .* dx

    dy = Base.circshift(u, [0, -1, 0]) - u
    dy[:,end,:] = dy[:,end-1,:]
    dy = ih[2] .* dy

    dz = Base.circshift(u, [0, 0, -1]) - u
    dz[:,:,end] = dz[:,:,end-1]
    dz = ih[3] .* dz;

    return [dx, dy, dz]
end


function validateinputs(u, mask, vsz)

    classes = {"single", "double"};
    attributes = {"real", "ndims", 3, "finite"};
    validateattributes(u, classes, attributes, mfilename, "u", 1);

    classes = {"logical", "numeric"};
    attributes = {"real", "ndims", 3, "size", size(u), "finite", "binary"};
    validateattributes(mask, classes, attributes, mfilename, "mask", 2);

    classes = {"numeric"};
    attributes = {"real", "vector", "numel", 3, "finite", ">", 0};
    validateattributes(vsz, classes, attributes, mfilename, "vsz", 3);

end
