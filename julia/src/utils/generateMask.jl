using Shell

include("../io/saveNii.jl")


function generateMask_(mag, vsz, betargs)

    magfile_temp = tempname(tempdir())
    magfile = [magfile_temp, ".nii"];
    magfile = join(magfile)
    saveNii(magfile, mag, vsz);

  # Generate the mask file using bet

    maskfile_temp = [magfile_temp, "_mask.nii.gz"]
    maskfile = join(maskfile_temp)

  # Add command to run bet in shell
    Shell.run("bet $magfile $magfile $betargs");
    Shell.run("gunzip $maskfile");
  # cmddd = ['bet ', magfile, ' ', magfile, ' ', betargs];

    maskfile_temp2 = [magfile_temp, "_mask.nii"]
    maskfile2 = join(maskfile_temp2)

    nii_obj = loadNii(maskfile2);

#   mask = logical(nii_obj.img);
    mask = nii_obj.img .!= 0

    return mask
end


function validateinputs(mag, vsz, betargs)

    if ndims(mag) < 4
        nd = 3;
    else
        nd = 4;
    end

    # classes = {"numeric"};
    # attributes = {"real", "ndims", nd};
    # validateattributes(mag, classes, attributes, mfilename, "mag", 1);

    # classes = {"numeric"};
    # attributes = {"real", "vector", "numel", 3, "finite", ">", 0};
    # validateattributes(vsz, classes, attributes, mfilename, "vsz", 2);

    # classes = {"char"};
    # validateattributes(betargs, classes, {}, mfilename, "betargs", 3);

end


function generateMask(mag, vsz, betargs)
#GENERATEMASK Automatic brain extraction using FSL's bet.
#
#   [mask] = GENERATEMASK(mag, vsz, [betargs]);
#
#   Inputs
#   ------
#       mag         magnitude image used as input for bet.
#       vsz         voxel size.
#
#       betargs     string containing command line arguments for bet.
#                   Default = '-m -n -f 0.5'.
#
#   Outputs
#   -------
#       mask        logical array containing binary mask.

    #<TODO> Find Julia equivalent narginchk
    #narginchk(2, 3);

    #<TODO> Pass these as default argiments from caller
    #if nargin < 3, betargs = '-m -n -f 0.5'; end

    # validateinputs(mag, vsz, betargs);

    mask = generateMask_(mag, vsz, betargs);

    return mask
end



