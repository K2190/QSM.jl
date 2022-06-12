# This will generate a mesh grid of dimesnions "dim" using the values
# in inp_array. This array also must be of size dim
function mesh_grid(xx, yy, zz)

  lenX = length(xx)
  lenY = length(yy)
  lenZ = length(zz)
  out_arrayX = zeros(Float64, lenY, lenX, lenZ)
  out_arrayY = zeros(Float64, lenY, lenX, lenZ)
  out_arrayZ = zeros(Float64, lenY, lenX, lenZ)

  for k in 1:length(zz)
    for i in 1:length(xx)
      for j in 1:length(yy)
        out_arrayX[j,i,k] = xx[i]
      end
    end
  end

  for k in 1:length(zz)
    for i in 1:length(xx)
      for j in 1:length(yy)
        out_arrayY[j,i,k] = yy[j]
      end
    end
  end

  for k in 1:length(zz)
    for i in 1:length(xx)
      for j in 1:length(yy)
        out_arrayZ[j,i,k] = zz[k]
      end
    end
  end

  return out_arrayX, out_arrayY, out_arrayZ
end



# temp = zeros(Float64, 2, 4, 2)
# N = size(temp)
# xx = range(-N[1]/2 , N[1]/2-1 , step=1) |> collect
# yy = range(-N[2]/2 , N[2]/2-1, step=1) |> collect
# zz = range(-N[3]/2 , N[3]/2-1, step=1) |> collect
#
# ky, kx, kz = mesh_grid(xx, yy, zz)
#
# println(kz)
