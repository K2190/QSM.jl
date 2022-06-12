 
using LinearAlgebra, Statistics

function jl_findNonZero(inp_array)
	cart_indicies = findall( x -> x != 0.0,  inp_array)
	lin_indices = LinearIndices(inp_array)
	non_zero_indicies = Vector{Int32}()
	for i = 1:length(cart_indicies)
	   index = lin_indices[cart_indicies[i]]
	   append!(non_zero_indicies, index)
	end
 
	return non_zero_indicies
end


function cropIndices(x)
#CROPINDICES [ix, iy, iz] = cropIndices(x)
#   See also CROP2MASK, UNCROP2MASK

    # <TODO> Check narginchk in Julia
    # narginchk(1, 1);

    s1 = sum(x, dims=1);
    s2 = sum(x, dims=2);

    tmparr1 = sum(s2, dims=3)
    ix = jl_findNonZero(tmparr1)

    tmparr2 = sum(s1, dims=3)
    iy = jl_findNonZero(tmparr2)

    tmparr3 = sum(s1, dims=2)
    iz = jl_findNonZero(tmparr3)

    ix = vec_jl(ix);
    iy = vec_jl(iy);
    iz = vec_jl(iz);

    return ix, iy, iz
end


# TEST cropIndices
#xx = [10.0 20.0 30.0]
#ix, iy, iz = cropIndices(xx)
#sz = [length(ix), length(iy), length(iz)]
#print(sz)


