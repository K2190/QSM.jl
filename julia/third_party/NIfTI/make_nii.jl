
include("../../src/utils/matlab_2_julia.jl")
include("../../src/utils/set_global_data.jl")

#  Make NIfTI structure specified by an N-D matrix. Usually, N is 3 for 
#  3D matrix [x y z], or 4 for 4D matrix with time series [x y z t]. 
#  Optional parameters can also be included, such as: voxel_size, 
#  origin, datatype, and description. 
#  
#  Once the NIfTI structure is made, it can be saved into NIfTI file 
#  using "save_nii" command (for more detail, type: help save_nii). 
#  
#  Usage: nii = make_nii(img, [voxel_size], [origin], [datatype], [description])
#
#  Where:
#
#	img:		Usually, img is a 3D matrix [x y z], or a 4D
#			matrix with time series [x y z t]. However,
#			NIfTI allows a maximum of 7D matrix. When the
#			image is in RGB format, make sure that the size
#			of 4th dimension is always 3 (i.e. [R G B]). In
#			that case, make sure that you must specify RGB
#			datatype, which is either 128 or 511.
#
#	voxel_size (optional):	Voxel size in millimeter for each
#				dimension. Default is [1 1 1].
#
#	origin (optional):	The AC origin. Default is [0 0 0].
#
#	datatype (optional):	Storage data type:
#		2 - uint8,  4 - int16,  8 - int32,  16 - float32,
#		32 - complex64,  64 - float64,  128 - RGB24,
#		256 - int8,  511 - RGB96,  512 - uint16,
#		768 - uint32,  1792 - complex128
#			Default will use the data type of 'img' matrix
#			For RGB image, you must specify it to either 128
#			or 511.
#
#	description (optional):	Description of data. Default is ''.
#
#  e.g.:
#     origin = [33 44 13]; datatype = 64;
#     nii = make_nii(img, [], origin, datatype);    % default voxel_size
#
#  NIFTI data format can be found on: http://nifti.nimh.nih.gov
#
#  - Jimmy Shen (jimmy@rotman-baycrest.on.ca)
#
function make_nii(img, vsz)

   nii_img = img  #varargin{1}
   dims = size(nii_img)

   # dims = size(nii.img);
   # dims = [length(dims) dims ones(1,8)];
   # dims = dims(1:8);

   tempDims = size(nii_img)
   dims = [tempDims[1] tempDims[2] tempDims[3]]
   dims = [length(dims) dims ones(1,8)];
   dims = convert(Array{Int16}, dims)
   dims = dims[1:8]

   # # dims = [length(dims) dims ones(1,8)]
   # # dims = dims[1:8]
   # new_dims = zeros(1)
   # append!(new_dims, length(dims))
   # for i=1:length(dims)
   #    append!(new_dims, dims[i])
   # end
   # for i=1:8
   #    append!(new_dims, 1)
   # end
   # dims = new_dims[1:8]
   

   voxel_size = [0, 1, 1, 1, 1, 1, 1, 1] #[0 ones(1,7)]
   voxel_size = convert_to_double(voxel_size)
   
   origin = zeros(1,5)
   descrip = ""

   # <TODO> Get a function to extract class
   nii_img_class = "single"

   if eltype(nii_img) == Float32
      nii_img_class = "single"
   elseif eltype(nii_img) == Float64
      nii_img_class = "double"
   else
      println("ERROR: Data type conversion is not valid")
   end

   # switch class(nii.img)
      if nii_img_class == "uint8"
         datatype = 2
      elseif nii_img_class == "int16"
         datatype = 4
      elseif nii_img_class == "int32"
         datatype = 8
      elseif nii_img_class == "single"
         if isreal(nii_img)
            datatype = 16
         else
            datatype = 32
         end
      elseif nii_img_class == "double"
         if isreal(nii_img)
            datatype = 64
         else
            datatype = 1792
         end
      elseif nii_img_class == "int8"
         datatype = 256
      elseif nii_img_class == "uint16"
         datatype = 512
      elseif nii_img_class == "uint32"
         datatype = 768
      else
         error("Datatype is not supported by make_nii.")
   end

   vsz_dbl = convert_to_double(vsz)
   if ~isempty(vsz)
         voxel_size[2] = vsz_dbl[1]
         voxel_size[3] = vsz_dbl[2]
         voxel_size[4] = vsz_dbl[3]
   end

   # if nargin > 1 & ~isempty(varargin{2})
   #    voxel_size(2:4) = convert_to_double(varargin{2})
   # end

   # if nargin > 2 & ~isempty(varargin{3})
   #    origin(1:3) = convert_to_double(varargin{3})
   # end

   # if nargin > 3 & ~isempty(varargin{4})
   #    datatype = convert_to_double(varargin{4})

   #    if datatype == 128 | datatype == 511
   #       dims(5) = []
   #       dims(1) = dims(1) - 1
   #       dims = [dims 1]
   #    end
   # end

   # if nargin > 4 & ~isempty(varargin{5})
   #    descrip = varargin{5}
   # end

   # if ndims(nii.img) > 7
   #    error("NIfTI only allows a maximum of 7 Dimension matrix.")
   # end

   temp1 = maximum(nii_img[:])
   maxval = round(convert_to_double(maximum(nii_img[:])))
   minval = round(convert_to_double(minimum(nii_img[:])))

   nii_hdr = make_header(dims, voxel_size, origin, datatype, descrip, maxval, minval);

   #<TODO> How to get the workspace info of nii#
   nii_hdr_dime_datatype = get_nii_obj_dime_datatype()

   # switch nii.hdr.dime.datatype
   if nii_hdr_dime_datatype == 2
      nii_img = convert_to_uint8(nii_img)
   elseif nii_hdr_dime_datatype == 4
      nii_img = convert_to_int16(nii_img)
   elseif nii_hdr_dime_datatype == 8
      nii_img = convert_to_int32(nii_img)
   elseif nii_hdr_dime_datatype == 16
      nii_img = convert_to_single(nii_img)
   elseif nii_hdr_dime_datatype == 32
      nii_img = convert_to_single(nii_img)
   elseif nii_hdr_dime_datatype == 64
      nii_img = convert_to_double(nii_img)
   elseif nii_hdr_dime_datatype == 128
      nii_img = convert_to_uint8(nii_img)
   elseif nii_hdr_dime_datatype == 256
      nii_img = convert_to_int8(nii_img)
   elseif nii_hdr_dime_datatype == 511
      img = convert_to_double(nii_img(:))
      img = convert_to_single((img - min(img))/(max(img) - min(img)))
      nii_img = reshape(img, size(nii_img))
      nii_hdr_dime_glmax = convert_to_double(max(img))
      nii_hdr_dime_glmin = convert_to_double(min(img))
   elseif nii_hdr_dime_datatype == 512
      nii_img = convert_to_uint16(nii_img)
   elseif nii_hdr_dime_datatype == 768
      nii_img = convert_to_uint32(nii_img)
   elseif nii_hdr_dime_datatype == 1792
      nii_img = convert_to_double(nii_img)
   else
      println("make_nii.jl: Datatype is not supported by make_nii.")
      ccall(:jl_exit, Cvoid, (Int32,), 86) 
   end

   # Definition of nii
# mutable struct nii
#    hdr::dsr
#    filetype::Int32	
#    fileprefix::String
#    img::Array{Float64}
# end

   nii_filetype = 2
   nii_fileprefix = ""
   nii_obj = nii(nii_hdr, nii_filetype, nii_fileprefix, nii_img)

   return nii_obj					# make_nii

end

# %---------------------------------------------------------------------
function make_header(dims, voxel_size, origin, datatype, descrip, maxval, minval)

   hdr_hk   = header_key_makenii()
   hdr_dime = image_dimension_makenii(dims, voxel_size, datatype, maxval, minval)
   hdr_hist = data_history_makenii(origin, descrip)
    
   dsr_obj = dsr(hdr_hk, hdr_dime, hdr_hist)
   return dsr_obj					# make_header
end

# %---------------------------------------------------------------------
function  header_key_makenii()

    hk_sizeof_hdr       = 348			# must be 348!
    hk_data_type        = codeunits("")
    hk_db_name          = codeunits("")
    hk_extents          = 0
    hk_session_error    = 0
    hk_regular          = 114   # ASCII value of "r" = 114
    hk_dim_info         = 0
    
    hk = header_key(hk_sizeof_hdr, hk_data_type, hk_db_name, hk_extents, hk_session_error, hk_regular, hk_dim_info)
    return hk				# header_key
end

# %---------------------------------------------------------------------
function image_dimension_makenii(dims, voxel_size, datatype, maxval, minval)
   
   dime_dim = dims
   dime_intent_p1 = 0
   dime_intent_p2 = 0
   dime_intent_p3 = 0
   dime_intent_code = 0
   dime_datatype = datatype

   # switch dime.datatype
   if dime_datatype == 2
      dime_bitpix = 8;  
      precision = "uint8"
   elseif dime_datatype == 4
      dime_bitpix = 16; 
      precision = "int16"
   elseif dime_datatype == 8
      dime_bitpix = 32; 
      precision = "int32"
   elseif dime_datatype == 16
      dime_bitpix = 32; 
      precision = "float32"
   elseif dime_datatype == 32
      dime_bitpix = 64; 
      precision = "float32"
   elseif dime_datatype == 64
      dime_bitpix = 64; 
      precision = "float64"
   elseif dime_datatype == 128
      dime_bitpix = 24;  
      precision = "uint8"
   elseif dime_datatype == 256 
      dime_bitpix = 8;  
      precision = "int8"
   elseif dime_datatype == 511
      dime_bitpix = 96;  
      precision = "float32"
   elseif dime_datatype == 512 
      dime_bitpix = 16; 
      precision = "uint16"
   elseif dime_datatype == 768 
      dime_bitpix = 32; 
      precision = "uint32"
   elseif dime_datatype == 1792
      dime_bitpix = 128; 
      precision = "float64"
   else
      error("Datatype is not supported by make_nii.")
   end
   
   dime_slice_start = 0
   dime_pixdim = voxel_size
   dime_vox_offset = 0
   dime_scl_slope = 0
   dime_scl_inter = 0
   dime_slice_end = 0
   dime_slice_code = 0
   dime_xyzt_units = 0
   dime_cal_max = 0
   dime_cal_min = 0
   dime_slice_duration = 0
   dime_toffset = 0
   dime_glmax = maxval
   dime_glmin = minval
   
   dime = image_dimension(dime_dim, dime_intent_p1, dime_intent_p2, dime_intent_p3, dime_intent_code,
                          dime_datatype, dime_bitpix, dime_slice_start, dime_pixdim, dime_vox_offset, dime_scl_slope,
                          dime_scl_inter, dime_slice_end, dime_slice_code, dime_xyzt_units, dime_cal_max,
                          dime_cal_min, dime_slice_duration, dime_toffset, dime_glmax, dime_glmin)

   return dime					# image_dimension
end

# %---------------------------------------------------------------------
function data_history_makenii(origin, descrip)
   
   uchar_data80 = zeros(UInt8, 80) 
   uchar_data24 = zeros(UInt8, 24) 
   uchar_data16 = zeros(UInt8, 16) 
   uchar_data4 = zeros(UInt8, 4) 

   hist_descrip = codeunits(descrip)
   hist_aux_file = codeunits("none")  #"none"
   hist_qform_code = 0
   hist_sform_code = 0
   hist_quatern_b = 0
   hist_quatern_c = 0
   hist_quatern_d = 0
   hist_qoffset_x = 0
   hist_qoffset_y = 0
   hist_qoffset_z = 0
   hist_srow_x = zeros(4)
   hist_srow_y = zeros(4)
   hist_srow_z = zeros(4)
   hist_intent_name = codeunits("")  #""
   hist_magic = codeunits("")  #""

   # hist_originator = zeros(3)  #convert(Array{Int16}, origin)
   # flip_orient = zeros(3)
   # rot_orient = zeros(3)

   hist = data_history(hist_descrip, hist_aux_file, hist_qform_code, hist_sform_code, hist_quatern_b, hist_quatern_c,
                       hist_quatern_d, hist_qoffset_x, hist_qoffset_y, hist_qoffset_z, hist_srow_x, hist_srow_y,
                       hist_srow_z, hist_intent_name, hist_magic)
                       #, hist_originator, flip_orient, rot_orient)
   
   return hist					# data_history
end
