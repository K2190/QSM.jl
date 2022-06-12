function dilateMask_(mask, vx)

    for ii = 1:vx
        mask = mask + bwperim(mask, 6)
    end

    return mask

end

function [mask] = dilateMask(mask, vx)
#DILATEMASK Dilate 3d binary mask.
#
#   [mask] = DILATEMASK(mask, vx);
#
#   Inputs
#   ------
#       mask    3d binary mask.
#       vx      number of voxels to dilate.
#
#   Outputs
#   -------
#       mask    dilated binary mask.
#
#   See also ERODEMASK

    #<TODO> Find Julia equivalent narginchk
    #narginchk(2, 2);

    if vx < 1
        return
    end

    #T = class(mask);
    T = typeof(mask[1])
    
    #mask = logical(mask);
    mask = !=(0).(mask)

    mask = dilateMask_(mask, vx);

    #mask = cast(mask, T);
    mask = convert(T, mask)

end



