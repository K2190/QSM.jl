
include("../utils/matlab_2_julia.jl")
include("../../third_party/NIfTI/load_nii.jl")

using Printf

function loadNii(filename)
#LOADNII Load NIfTI dataset.
#   Wrapper for Jimmy Shen's NIfTI and ANALYZE toolbox.
#
#   [nii] = LOADNII(filename);
#
#   See also LOAD_NII, MAKE_NII, SAVE_NII, SAVENII

    #narginchk(1, 1);

    println(@sprintf("Loading file %s", filename))

    if jl_exist(filename, "file") == 2
        _, _, ext = jl_fileparts(filename)

        if ext == ".gz"
            gunzip(filename)
            nii = load_nii(filename[1:end-3]);
            delete(filename[1:end-3])
            return nii
        elseif ext == ".nii"
            nii = load_nii(filename)
            return nii
        else
            println("wrong file extension %s", filename)
        end
    else
        println(@sprintf("file %s does not exist", filename))
    end

    return nii
end
