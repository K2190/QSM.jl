function crop2mask(x, mask)
#CROP2MASK [u, ix, iy, iz] = crop2mask(x, [mask])
#   See also UNCROP2MASK, CROPINDICES

    #<TODO> Find Julia equivalent narginchk
    #narginchk(1, 2)
    #if nargin < 2, mask = x; end

    ix, iy, iz = cropIndices(mask)
    u = x[ix,iy,iz]

    return [u, ix, iy, iz] 
end
