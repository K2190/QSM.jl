function [dx, dy, dz] = gradc(u, mask, vsz)
#GRADC First order centered difference gradient.
#   Boundary condition: u^{k+1} - 2u^{k} + u^{k-1} = 0.
#
#   [dx, dy, dz] = GRADC(u, [mask], [vsz]);
#
#   See also GRADF, GRADB

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
        gradc_mex(dx, dy, dz, u, vsz);
    else
        gradcm_mex(dx, dy, dz, u, mask, vsz);
    end

    return [dx, dy, dz]
end


function grad_(u, h)

    ih = 1 ./ h;

    dx = Base.circshift(u, [-1,0,0]) - Base.circshift(u, [1,0,0]);
    dx = 0.5 .* ih[1] .* dx;
    dx[1,:,:] = ih[1] .* (u[2,:,:] - u[1,:,:]);
    dx[end,:,:] = ih[1] .* (u[end,:,:] - u[end-1,:,:]);

    dy = Base.circshift(u, [0,-1,0]) - Base.circshift(u, [0,1,0]);
    dy = 0.5 .* ih[2] .* dy;
    dy[:,1,:] = ih[2] .* (u[:,2,:] - u[:,1,:]);
    dy[:,end,:] = ih[2] .* (u[:,end,:] - u[:,end-1,:]);

    dz = Base.circshift(u, [0,0,-1]) - Base.circshift(u, [0,0,1]);
    dz = 0.5 .* ih[3] .* dz;
    dz[:,:,1] = ih[3] .* (u[:,:,2] - u[:,:,1]);
    dz[:,:,end] = ih[3] .* (u[:,:,end] - u[:,:,end-1]);

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
