
include("load_nii_hdr.jl")
include("load_nii_img.jl")
include("xform_nii.jl")

mutable struct nii
    hdr::dsr
	 filetype::Int32	
	 fileprefix::String
    img::Array{Float64}
end


#  Load NIFTI or ANALYZE dataset. Support both *.nii and *.hdr/*.img
#  file extension. If file extension is not provided, *.hdr/*.img will
#  be used as default.
#
#  A subset of NIFTI transform is included. For non-orthogonal rotation,
#  shearing etc., please use 'reslice_nii.m' to reslice the NIFTI file.
#  It will not cause negative effect, as long as you remember not to do
#  slice time correction after reslicing the NIFTI file. Output variable
#  nii will be in RAS orientation, i.e. X axis from Left to Right,
#  Y axis from Posterior to Anterior, and Z axis from Inferior to
#  Superior.
#  
#  Usage: nii = load_nii(filename, [img_idx], [dim5_idx], [dim6_idx], ...
#			[dim7_idx], [old_RGB], [tolerance], [preferredForm])
#  
#  filename  - 	NIFTI or ANALYZE file name.
#  
#  img_idx (optional)  -  a numerical array of 4th dimension indices,
#	which is the indices of image scan volume. The number of images
#	scan volumes can be obtained from get_nii_frame.m, or simply
#	hdr.dime.dim(5). Only the specified volumes will be loaded. 
#	All available image volumes will be loaded, if it is default or
#	empty.
#
#  dim5_idx (optional)  -  a numerical array of 5th dimension indices.
#	Only the specified range will be loaded. All available range
#	will be loaded, if it is default or empty.
#
#  dim6_idx (optional)  -  a numerical array of 6th dimension indices.
#	Only the specified range will be loaded. All available range
#	will be loaded, if it is default or empty.
#
#  dim7_idx (optional)  -  a numerical array of 7th dimension indices.
#	Only the specified range will be loaded. All available range
#	will be loaded, if it is default or empty.
#
#  old_RGB (optional)  -  a scale number to tell difference of new RGB24
#	from old RGB24. New RGB24 uses RGB triple sequentially for each
#	voxel, like [R1 G1 B1 R2 G2 B2 ...]. Analyze 6.0 from AnalyzeDirect
#	uses old RGB24, in a way like [R1 R2 ... G1 G2 ... B1 B2 ...] for
#	each slices. If the image that you view is garbled, try to set 
#	old_RGB variable to 1 and try again, because it could be in
#	old RGB24. It will be set to 0, if it is default or empty.
#
#  tolerance (optional) - distortion allowed in the loaded image for any
#	non-orthogonal rotation or shearing of NIfTI affine matrix. If 
#	you set 'tolerance' to 0, it means that you do not allow any 
#	distortion. If you set 'tolerance' to 1, it means that you do 
#	not care any distortion. The image will fail to be loaded if it
#	can not be tolerated. The tolerance will be set to 0.1 (10%), if
#	it is default or empty.
#
#  preferredForm (optional)  -  selects which transformation from voxels
#	to RAS coordinates; values are s,q,S,Q.  Lower case s,q indicate
#	"prefer sform or qform, but use others if preferred not present". 
#	Upper case indicate the program is forced to use the specificied
#	tranform or fail loading.  'preferredForm' will be 's', if it is
#	default or empty.	- Jeff Gunter
#
#  Returned values:
#  
#  nii structure:
#
#	hdr -		struct with NIFTI header fields.
#
#	filetype -	Analyze format .hdr/.img (0); 
#			NIFTI .hdr/.img (1);
#			NIFTI .nii (2)
#
#	fileprefix - 	NIFTI filename without extension.
#
#	machine - 	machine string variable.
#
#	img - 		3D (or 4D) matrix of NIFTI data.
#
#	original -	the original header before any affine transform.
#  
#  Part of this file is copied and modified from:
#  http://www.mathworks.com/matlabcentral/fileexchange/1878-mri-analyze-tools
#  
#  NIFTI data format can be found on: http://nifti.nimh.nih.gov
#  
#  - Jimmy Shen (jimmy@rotman-baycrest.on.ca)
#
function load_nii(filename, img_idx=[], dim5_idx=[], dim6_idx=[], dim7_idx=[], 
                  old_RGB=0.0, tolerance=0.1, preferredForm='s')

   # if ~exist("filename","var")
   #    error("Usage: nii = load_nii(filename, [img_idx], [dim5_idx], [dim6_idx], [dim7_idx], [old_RGB], [tolerance], [preferredForm])");
   # end

   if isempty(img_idx)
      img_idx = []
   end

   if isempty(dim5_idx)
      dim5_idx = []
   end

   if isempty(dim6_idx)
      dim6_idx = []
   end

   if isempty(dim7_idx)
      dim7_idx = []
   end

   if isempty(old_RGB)
      old_RGB = 0
   end

   if isempty(tolerance)
      tolerance = 0.1			# 10 percent
   end

   if isempty(preferredForm)
      preferredForm= "s"		# Jeff
   end

   # Not needed in Julia, just hard code it
   #v = version
   v= 7.8

   #  Check file extension. If .gz, unpack it into temp folder
   #

   #f_name = filename[end-2:end]

   _, _, extn = jl_fileparts(filename)

   if extn != ".nii"  # Temporary, remove later

   if length(filename) > 2 && filename[end-2:end] == ".gz"
      if filename[end-6:end] != ".img.gz" && filename[end-6:end] != ".hdr.gz" && filename[end-6:end] != ".nii.gz"
         println("Error: Please check filename.")
      end

      
      if parse(Float64, v[1:3]) < 7.1 
         println("Please use MATLAB 7.1 (with java) and above, or run gunzip outside MATLAB.")
      elseif filename[end-6:end] == ".img.gz"
         filename1 = filename
         filename2 = filename
         filename2[end-6:end] = ""
         filename2 = [filename2, ".hdr.gz"]

         tmpDir = tempname
         mkdir(tmpDir)
         gzFileName = filename

         filename1 = gunzip(filename1, tmpDir)
         filename2 = gunzip(filename2, tmpDir)
         filename = char(filename1)	# convert from cell to string
      elseif filename[end-6:end] == ".hdr.gz"
         filename1 = filename
         filename2 = filename
         filename2[end-6:end] = ""
         filename2 = [filename2, ".img.gz"]

         tmpDir = tempname
         mkdir(tmpDir)
         gzFileName = filename

         filename1 = gunzip(filename1, tmpDir)
         filename2 = gunzip(filename2, tmpDir)
         filename = char(filename1)	# convert from cell to string
      elseif filename[end-6:end] == ".nii.gz"
         tmpDir = tempname
         mkdir(tmpDir)
         gzFileName = filename
         filename = gunzip(filename, tmpDir)
         filename = char(filename)	# convert from cell to string
      end
   end

   end   # Temporary, remove later
   #  Read the dataset header
   #
   
   nii_hdr, nii_filetype, nii_fileprefix = load_nii_hdr(filename)

   hdr_hist_srow_x = nii_hdr.hist.srow_x
   hdr_hist_srow_y = nii_hdr.hist.srow_y
   hdr_hist_srow_z = nii_hdr.hist.srow_z

   #  Read the header extension
   #
   # nii.ext = load_nii_ext(filename);

   #  Read the dataset body
   nii_obj_img, nii_hdr_2 = load_nii_img(nii_hdr,  nii_filetype,  nii_fileprefix,
		                                     "ieee",img_idx,dim5_idx,dim6_idx,dim7_idx,old_RGB)

   nii_obj = nii(nii_hdr, nii_filetype, nii_fileprefix, nii_obj_img)

   # #  Perform some of sform/qform transform
   # #

   hdr_hist_srow_x = nii_obj.hdr.hist.srow_x
   hdr_hist_srow_y = nii_obj.hdr.hist.srow_y
   hdr_hist_srow_z = nii_obj.hdr.hist.srow_z

   nii_obj2 = xform_nii(nii_obj, tolerance, preferredForm)

   # mag = nii_obj2.img[:,:,:,:,1]

   # #  Clean up after gunzip
   # #
   # if exist("gzFileName", "var")

   #    #  fix fileprefix so it doesn't point to temp location
   #    #
   #    nii.fileprefix = gzFileName[1:end-7]
   #    rmdir(tmpDir,"s")
   # end

   return nii_obj
end					# load_nii

