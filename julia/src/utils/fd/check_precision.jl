function check_precision(var, type)

	if isa(var, Array)
		if lowercase(type) == "single"
			if typeof(var[1]) == Int32 || typeof(var[1]) == Float32
				return true;
			end
		end

		if lowercase(type) == "double"
			if typeof(var[1]) == Int64 || typeof(var[1]) == Float64
				return true;
			end
		end

	else
		if lowercase(type) == "single"
			if typeof(var) == Int32 || typeof(var) == Float32
				return true;
			end
		end

		if lowercase(type) == "double"
			if typeof(var) == Int64 || typeof(var) == Float64
				return true;
			end
		end
	end

	return false
end

# xx = [ [1; 2] [3; 4]]
# println(check_precision(xx, "single"))

# dx = convert(Matrix{Float32}, xx)
# println(check_precision(dx, "single"))

# println(check_precision(5, "single"))

# println(typeof(5))

# println(isinteger(5.2))