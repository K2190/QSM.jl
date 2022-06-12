# This will generate a mesh grid of dimesnions "dim" using the values
# in inp_array. This array also must be of size dim
# Equivalent to ndgrid of MATLAB
function nd_grid(xx, yy, zz)

  lenX = length(xx)
  lenY = length(yy)
  lenZ = length(zz)
  out_arrayX = zeros(Float64, lenX, lenY, lenZ)
  out_arrayY = zeros(Float64, lenX, lenY, lenZ)
  out_arrayZ = zeros(Float64, lenX, lenY, lenZ)

  out_arrayX2 = zeros(Float64, lenX, lenY)
  out_arrayY2 = zeros(Float64, lenX, lenY)

  for i in 1:length(xx)
    for j in 1:length(yy)
      out_arrayX2[i,j] = xx[i]
    end
  end
  for k in 1:length(zz)
      out_arrayX[:, :, k] = out_arrayX2
  end


  for i in 1:length(xx)
    for j in 1:length(yy)
      out_arrayY2[i,j] = yy[j]
    end
  end
  for k in 1:length(zz)
    out_arrayY[:, :, k] = out_arrayY2
  end


  for i in 1:length(xx)
    for j in 1:length(yy)
      for k in 1:length(zz)
        out_arrayZ[i,j,k] = zz[k]
      end
    end
  end

  return out_arrayX, out_arrayY, out_arrayZ
end



# TEST CODE
# start = -2
# end1 = 1
# x = collect(start : 1 : end1)
# y = collect(start : 1 : end1)
# z = collect(start : 1 : end1)

# x = collect(start : 1 : end1)
# y = collect(start-1 : 1 : end1)
# z = collect(start-2 : 1 : end1)

# X,Y,Z = mesh_grid2(x,y,z)

# println((x))
# println((X))
# println((Y))
# println((Z))
