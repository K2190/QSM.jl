
using ZipFile

#================================================================================#
function zeros_like(sz, str, inp_var)
    out_var = zeros(sz)
    if isa(inp_var, Array)
       out_var = convert(Array{typeof(inp_var[1])}, out_var)
    else
       out_var = convert(Array{typeof(inp_var)}, out_var)
    end
    return out_var
 end

 #================================================================================#
function convert_to_typeT(in_var, T)
    if isa(in_var, Array)
        out_var = convert(Array{T}, in_var)
    elseif isa(in_var, Matrix)
        out_var = convert(Matrix{T}, in_var)
    elseif isa(in_var, Vector)
        out_var = convert(Vector{T}, in_var)
    else
        out_var = convert(T, in_var)
    end
    
    return out_var
end

#================================================================================#
function convert_to_uint8(in_var)
    if isa(in_var, Array)
        out_var = convert(Array{UInt8}, in_var)
    elseif isa(in_var, Matrix)
        out_var = convert(Matrix{UInt8}, in_var)
    elseif isa(in_var, Vector)
        out_var = convert(Vector{UInt8}, in_var)
    else
        out_var = convert(UInt8, in_var)
    end
    
    return out_var
end

#================================================================================#
function convert_to_double(in_var)
    # Matlab equivalent of double
    # https://in.mathworks.com/help/symbolic/double.html
    
        if isa(in_var, Array)
            out_var = convert(Array{Float64}, in_var)
        elseif isa(in_var, Matrix)
            out_var = convert(Matrix{Float64}, in_var)
        elseif isa(in_var, Vector)
            out_var = convert(Vector{Float64}, in_var)
        else
            out_var = convert(Float64, in_var)
        end
    
    return out_var
end

#================================================================================#
function convert_to_single(in_var)
    # Matlab equivalent of double
    # https://in.mathworks.com/help/symbolic/double.html
    
        if isa(in_var, Array)
            out_var = convert(Array{Float32}, in_var)
        elseif isa(in_var, Matrix)
            out_var = convert(Matrix{Float32}, in_var)
        elseif isa(in_var, Vector)
            out_var = convert(Vector{Float32}, in_var)
        else
            out_var = convert(Float32, in_var)
        end
    
    return out_var
end

#================================================================================#
function jl_exist(variable, type)
# Equivalent of exist function in matlab
# https://in.mathworks.com/help/matlab/ref/exist.html

    if type == "file"
        if isfile(variable)
            return 2
        end
    end

    return -1
end

#================================================================================#
function jl_fileparts(file)
# Equivalent of fileparts function in matlab
# https://in.mathworks.com/help/matlab/ref/fileparts.html
    filepath = getFileDirectory(file)
    name  = getFileName(file)
    ext = getFileExtension(file)

    return [filepath, name, ext]
end

#================================================================================#
function getFileExtension(filename)
    fileExtn = filename[findlast(isequal('.'),filename):end]
    return fileExtn
end

#================================================================================#
function getFileName(filename)
    return basename(filename)
end

#================================================================================#
function getFileDirectory(filename)
    return dirname(filename)
end

#================================================================================#

function jl_unzip(file,exdir="")
    fileFullPath = isabspath(file) ?  file : joinpath(pwd(),file)
    basePath = dirname(fileFullPath)
    outPath = (exdir == "" ? basePath : (isabspath(exdir) ? exdir : joinpath(pwd(),exdir)))
    isdir(outPath) ? "" : mkdir(outPath)
    zarchive = ZipFile.Reader(fileFullPath)
    for f in zarchive.files
        fullFilePath = joinpath(outPath,f.name)
        if (endswith(f.name,"/") || endswith(f.name,"\\"))
            mkdir(fullFilePath)
        else
            write(fullFilePath, read(f))
        end
    end
    close(zarchive)
end

#================================================================================#
####  TESTING CODE

# println(getFileExtension(fname1))
# println(getFileName(fname1))
# println(getFileDirectory(fname1))

# jl_unzip(fname1, "mydata")
# zarchive = ZipFile.Reader(fname1)
# myf = zarchive.files[1]
# dir1 = getFileDirectory(fname1)
# parts = [dir1, myf.name]
# new_file = joinpath(parts...)
# out =  open(new_file,"w")
# println(new_file)
# write(out,myf) # error
# close(out)
# close(zarchive)



