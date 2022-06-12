include("../utils/matlab_2_julia.jl")
include("../../third_party/NIfTI/make_nii.jl")
include("../../third_party/NIfTI/save_nii.jl")

function saveNii(filename, img, vsz)
#SAVENII Save NIfTI dataset.
#   Wrapper for Jimmy Shen's NIfTI and ANALYZE toolbox.
#
#   [nii] = SAVENII(filename, img, [vsz]);
#
#   See also MAKE_NII, SAVE_NII, LOAD_NII, LOADNII

    #narginchk(2, 3);

    #if nargin < 3, vsz = [1, 1, 1]; end


    pathstr, _, _ = jl_fileparts(filename)

    
    if ~isdir(pathstr) #exist(pathstr, "dir") != 7
        status, msg = mkdir(pathstr)
        if status != 1
            println(msg)
        end
    end

    nii_obj = make_nii(img, vsz)
    save_nii(nii_obj, filename)

    return nii_obj
end
