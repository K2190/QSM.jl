
using AbstractFFTs

include("../mesh_grid.jl")
include("../matlab_2_julia.jl")

function dipoleKernel(sz, vsz, bdir, space, T)
#DIPOLEKERNEL Generate kspace dipole kernel.
#   The dipole kernel can be constructed in image space or in k-space. The
#   kspace kernel will be shifted such that the DC component is at index
#   (1,1,1).
#
#   [D] = DIPOLEKERNEL(sz, vsz, [bdir], [space], [T]);
#
#   Inputs
#   ------
#       sz      image size (3 element vector).
#       vsz     voxel size (3 element vector).
#
#       bdir    unit vector of B field direction.
#               default: [0, 0, 1]
#       space   'i'magespace or 'k'space.
#               default: 'k'
#       T       floating point precision. 'single' or 'double'.
#               default: 'double'
#
#   Outputs
#   -------
#       D       dipole kernel.

    # narginchk(2, 5);

    # if nargin < 5, T = "double"; end
    # if nargin < 4 || isempty(space), space = "k"; end
    # if nargin < 3 || isempty(bdir), bdir = [0, 0, 1]; end

    scalar_val = ~(isa(vsz, Array) || isa(vsz, Vector) || isa(vsz, Matrix))
    if scalar_val
        vsz = [vsz, vsz, vsz]
    end

    # if norm(bdir) ~= 1
    #     warning("norm(bdir) = %g ~= 1, normalizing...")
    #     bdir = bdir ./ norm(bdir)
    # end

    sz = reshape(collect(sz), 1, :)
    vsz = reshape(vsz, 1, :)

    # switch lower(space(1))
    if space[1] == "i"
            X, Y, Z = dipoleGrid(sz, vsz)

            rz = bdir(1).*X + bdir(2).*Y + bdir(3).*Z
            r2  = X.*X + Y.*Y + Z.*Z;

            den = 4*pi .* sqrt(r2 .* r2 .* r2 .* r2 .* r2)

            D = (3*(rz.*rz) - r2) ./ den
            D[den == 0] = 0
    else
            temp = sz.*vsz
            temp2 = zeros(eltype(temp), size(temp))
            temp2 = 1 ./temp
            X, Y, Z = dipoleGrid(sz, 1 ./(sz.*vsz))

            kz = bdir[1] .* X + bdir[2] .* Y + bdir[3] .* Z
            k2 = X.*X + Y.*Y + Z.*Z

            D = 1/3 .- (kz.*kz ./ k2)
            D1 = (kz.*kz ./ k2)

            D = ifftshift(D)
            D[1] = 0

    end

    D = convert_to_typeT(D, T)
    return D
end


function dipoleGrid(sz, h)

    # Making sure both even and odd sizes will be centered properly
    # low  = floor(sz / 2);
    # high = low - ~mod(sz, 2);

    low  = floor.(sz./ 2)
    msz = mod.(sz, 2)
    msz2 = replace(msz, 0 =>1, 1 =>0)
    high = low .- msz2

    low  = h .* low
    high = h .* high

    xx = range(-low[1] , high[1], step=h[1]) |> collect
    yy = range(-low[2] , high[2], step=h[2]) |> collect
    zz = range(-low[3] , high[3], step=h[3]) |> collect
    X, Y, Z = mesh_grid(xx, yy, zz)

    # X, Y, Z = ndgrid(-low[1] : h[1] : high[1],
    #                    -low[2] : h[2] : high[2],
    #                    -low[3] : h[3] : high[3]);

    return [X, Y, Z]
end
