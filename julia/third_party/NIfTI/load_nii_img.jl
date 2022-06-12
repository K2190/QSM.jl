#  internal function

#  - Jimmy Shen (jimmy@rotman-baycrest.on.ca)
function jl_sub2ind(inp_sz, n1, n2, n3, n4, n5, n6, n7)
    dummy_arr = zeros(Int8, inp_sz...)
    index = dummy_arr[CartesianIndex(n1, n2, n3, n4, n5, n6, n7)]
    return index
end

function load_nii_img(hdr,filetype,fileprefix,machine,img_idx,dim5_idx,dim6_idx,dim7_idx,old_RGB)

   # if ~exist("hdr","var") | ~exist("filetype","var") | ~exist("fileprefix","var") | ~exist("machine","var")
   #    error("Usage: [img,hdr] = load_nii_img(hdr,filetype,fileprefix,machine,[img_idx],[dim5_idx],[dim6_idx],[dim7_idx],[old_RGB]);");
   # end

   if isempty(img_idx) | hdr.dime.dim[5]<1
      img_idx = []
   end

   if  isempty(dim5_idx) | hdr.dime.dim[6]<1
      dim5_idx = []
   end

   if isempty(dim6_idx) | hdr.dime.dim[7]<1
      dim6_idx = []
   end

   if isempty(dim7_idx) | hdr.dime.dim[8]<1
      dim7_idx = []
   end

   if isempty(old_RGB)
      old_RGB = 0
   end

   #  check img_idx
   #
   if !isempty(img_idx) && ~isnumeric(img_idx)
      println("img_idx should be a numerical array.")
   end

   if length(unique(img_idx)) != length(img_idx)
      println("Duplicate image index in " + img_idx)
   end

   if ~isempty(img_idx) && (min(img_idx) < 1 || max(img_idx) > hdr.dime.dim[5])
      max_range = hdr.dime.dim[5]

      if max_range == 1
         println([img_idx + " should be 1."])
      else
         range = ["1 " num2str(max_range)]
         println([img_idx + " should be an integer within the range of [" range "]."])
      end
   end

   #  check dim5_idx
   #
   if ~isempty(dim5_idx) && ~isnumeric(dim5_idx)
      println(dim5_idx + " should be a numerical array.")
   end

   if length(unique(dim5_idx)) != length(dim5_idx)
      println("Duplicate index in " + dim5_idx)
   end

   if ~isempty(dim5_idx) && (min(dim5_idx) < 1 || max(dim5_idx) > hdr.dime.dim[6])
      max_range = hdr.dime.dim[6]

      if max_range == 1
         println([dim5_idx + " should be 1."])
      else
         range = ["1 " num2str(max_range)]
         println([dim5_idx + " should be an integer within the range of [" range "]."])
      end
   end

   #  check dim6_idx
   #
   if ~isempty(dim6_idx) && ~isnumeric(dim6_idx)
      println(dim6_idx + " should be a numerical array.")
   end

   if length(unique(dim6_idx)) != length(dim6_idx)
      println("Duplicate index in " + dim6_idx)
   end

   if ~isempty(dim6_idx) && (min(dim6_idx) < 1 || max(dim6_idx) > hdr.dime.dim[7])
      max_range = hdr.dime.dim[7]

      if max_range == 1
         println([dim6_idx + " should be 1."])
      else
         range = ["1 " num2str(max_range)]
         println([dim6_idx + " should be an integer within the range of [" range "]."])
      end
   end

   #  check dim7_idx
   #
   if ~isempty(dim7_idx) && ~isnumeric(dim7_idx)
      println("dim7_idx should be a numerical array.")
   end

   if length(unique(dim7_idx)) != length(dim7_idx)
      println("Duplicate index in " + dim7_idx)
   end

   if ~isempty(dim7_idx) && (min(dim7_idx) < 1 || max(dim7_idx) > hdr.dime.dim(8))
      max_range = hdr.dime.dim(8)

      if max_range == 1
         println([dim7_idx + " should be 1."])
      else
         range = ["1 " num2str(max_range)]
         println([dim7_idx + " should be an integer within the range of [" range "]."])
      end
   end

   img, hdr = read_image(hdr,filetype,fileprefix,machine,img_idx,dim5_idx,dim6_idx,dim7_idx,old_RGB)

   return	[img, hdr]			# load_nii_img
end

function read_image(hdr,filetype,fileprefix,machine,img_idx,dim5_idx,dim6_idx,dim7_idx,old_RGB)

   if filetype == 0 || filetype == 1
      fn = [fileprefix ".img"]
   end
   if filetype == 2
      fn = [fileprefix ".nii"]
   end

   filename = join(fn)
   # fid = fopen(fn,"r",machine)

   # Open the IO stream of the file
   open(filename) do io

   # if fid < 0
   #    msg = sprintf("Cannot open file %s.",fn)
   #    error(msg)
   # end

   #  Set bitpix according to datatype
   #
   #  /*Acceptable values for datatype are*/ 
   #
   #     0 None                     (Unknown bit per voxel) # DT_NONE, DT_UNKNOWN 
   #     1 Binary                         (ubit1, bitpix=1) # DT_BINARY 
   #     2 Unsigned char         (uchar or uint8, bitpix=8) # DT_UINT8, NIFTI_TYPE_UINT8 
   #     4 Signed short                  (int16, bitpix=16) # DT_INT16, NIFTI_TYPE_INT16 
   #     8 Signed integer                (int32, bitpix=32) # DT_INT32, NIFTI_TYPE_INT32 
   #    16 Floating point    (single or float32, bitpix=32) # DT_FLOAT32, NIFTI_TYPE_FLOAT32 
   #    32 Complex, 2 float32      (Use float32, bitpix=64) # DT_COMPLEX64, NIFTI_TYPE_COMPLEX64
   #    64 Double precision  (double or float64, bitpix=64) # DT_FLOAT64, NIFTI_TYPE_FLOAT64 
   #   128 uint8 RGB                 (Use uint8, bitpix=24) # DT_RGB24, NIFTI_TYPE_RGB24 
   #   256 Signed char            (schar or int8, bitpix=8) # DT_INT8, NIFTI_TYPE_INT8 
   #   511 Single RGB              (Use float32, bitpix=96) # DT_RGB96, NIFTI_TYPE_RGB96
   #   512 Unsigned short               (uint16, bitpix=16) # DT_UNINT16, NIFTI_TYPE_UNINT16 
   #   768 Unsigned integer             (uint32, bitpix=32) # DT_UNINT32, NIFTI_TYPE_UNINT32 
   #  1024 Signed long long              (int64, bitpix=64) # DT_INT64, NIFTI_TYPE_INT64
   #  1280 Unsigned long long           (uint64, bitpix=64) # DT_UINT64, NIFTI_TYPE_UINT64 
   #  1536 Long double, float128  (Unsupported, bitpix=128) # DT_FLOAT128, NIFTI_TYPE_FLOAT128 
   #  1792 Complex128, 2 float64  (Use float64, bitpix=128) # DT_COMPLEX128, NIFTI_TYPE_COMPLEX128 
   #  2048 Complex256, 2 float128 (Unsupported, bitpix=256) # DT_COMPLEX128, NIFTI_TYPE_COMPLEX128 
   #

   if hdr.dime.datatype == 1
      hdr.dime.bitpix = 1
      precision = "ubit1"
   elseif hdr.dime.datatype == 2
      hdr.dime.bitpix = 8
      precision = "uint8"
   elseif hdr.dime.datatype == 4
      hdr.dime.bitpix = 16
      precision = "int16"
   elseif hdr.dime.datatype == 8
      hdr.dime.bitpix = 32
      precision = "int32"
   elseif hdr.dime.datatype == 16
      hdr.dime.bitpix = 32
      precision = "float32"
   elseif hdr.dime.datatype == 32
      hdr.dime.bitpix = 64
      precision = "float32"
   elseif hdr.dime.datatype == 64
      hdr.dime.bitpix = 64
      precision = "float64"
   elseif hdr.dime.datatype == 128
      hdr.dime.bitpix = 24
      precision = "uint8"
   elseif hdr.dime.datatype == 256
      hdr.dime.bitpix = 8
      precision = "int8"
   elseif hdr.dime.datatype == 511
      hdr.dime.bitpix = 96
      precision = "float32"
   elseif hdr.dime.datatype == 512
      hdr.dime.bitpix = 16
      precision = "uint16"
   elseif hdr.dime.datatype == 768
      hdr.dime.bitpix = 32
      precision = "uint32"
   elseif hdr.dime.datatype == 1024
      hdr.dime.bitpix = 64
      precision = "int64"
   elseif hdr.dime.datatype == 1280
      hdr.dime.bitpix = 64
      precision = "uint64"
   elseif hdr.dime.datatype == 1792
      hdr.dime.bitpix = 128
      precision = "float64"
   else
      println("This datatype is not supported") 
      ccall(:jl_exit, Cvoid, (Int32,), 86) 
   end

   # Check this LATER
   # hdr.dime.dim(find(hdr.dime.dim < 1)) = 1

   #  move pointer to the start of image block
   #
   if filetype == 0 || filetype == 1
      seek(io, 0)
   elseif filetype == 2
      seek(io, Int(hdr.dime.vox_offset))
   end

   #  Load whole image block for old Analyze format or binary image;
   #  otherwise, load images that are specified in img_idx, dim5_idx,
   #  dim6_idx, and dim7_idx
   #
   #  For binary image, we have to read all because pos can not be
   #  seeked in bit and can not be calculated the way below.
   #

   arr_ones = [1, 1, 1, 1]
   if hdr.dime.datatype == 1 || hdr.dime.dim[5:8] == arr_ones || 
	   (isempty(img_idx) && isempty(dim5_idx) && isempty(dim6_idx) && isempty(dim7_idx))

      #  For each frame, precision of value will be read 
      #  in img_siz times, where img_siz is only the 
      #  dimension size of an image, not the byte storage
      #  size of an image.
      #

      img_siz = prod(hdr.dime.dim[2:8])

      #  For complex float32 or complex float64, voxel values
      #  include [real, imag]
      #
      if hdr.dime.datatype == 32 | hdr.dime.datatype == 1792
         img_siz = img_siz * 2
      end
	 
      #MPH: For RGB24, voxel values include 3 separate color planes
      #
      if hdr.dime.datatype == 128 | hdr.dime.datatype == 511
	      img_siz = img_siz * 3
      end

      println("Reading image data")
      
      img = zeros(Float32, img_siz)

      if precision == "int16"
         img = convert(Array{Int16}, img)
         println("Image data format used: INT16")
      elseif precision == "float64"
         img = convert(Array{Float64}, img)
         println("Image data format used: Float64")
      elseif precision == "float32"
         img = convert(Array{Float32}, img)
         println("Image data format used: Float32")
      else
         println("ERROR: Data type conversion is not valid")
         ccall(:jl_exit, Cvoid, (Int32,), 86) 
         return
      end

      #read(io, img, sprintf("*%s",precision))
      read!(io, img)

      d1 = hdr.dime.dim[2]
      d2 = hdr.dime.dim[3]
      d3 = hdr.dime.dim[4]
      d4 = hdr.dime.dim[5]
      d5 = hdr.dime.dim[6]
      d6 = hdr.dime.dim[7]
      d7 = hdr.dime.dim[8]

      if isempty(img_idx)
         img_idx = 1:d4
      end

      if isempty(dim5_idx)
         dim5_idx = 1:d5
      end

      if isempty(dim6_idx)
         dim6_idx = 1:d6
      end

      if isempty(dim7_idx)
         dim7_idx = 1:d7
      end
   else

      d1 = hdr.dime.dim[2]
      d2 = hdr.dime.dim[3]
      d3 = hdr.dime.dim[4]
      d4 = hdr.dime.dim[5]
      d5 = hdr.dime.dim[6]
      d6 = hdr.dime.dim[7]
      d7 = hdr.dime.dim[8]

      if isempty(img_idx)
         img_idx = 1:d4
      end

      if isempty(dim5_idx)
         dim5_idx = 1:d5
      end

      if isempty(dim6_idx)
         dim6_idx = 1:d6
      end

      if isempty(dim7_idx)
         dim7_idx = 1:d7
      end

      #  compute size of one image
      #
      img_siz = prod(hdr.dime.dim[2:4])

      #  For complex float32 or complex float64, voxel values
      #  include [real, imag]
      #
      if hdr.dime.datatype == 32 | hdr.dime.datatype == 1792
         img_siz = img_siz * 2
      end

      #MPH: For RGB24, voxel values include 3 separate color planes
      #
      if hdr.dime.datatype == 128 | hdr.dime.datatype == 511
         img_siz = img_siz * 3
      end

      # preallocate img
      img = zeros(img_siz, length(img_idx)*length(dim5_idx)*length(dim6_idx)*length(dim7_idx) )
      currentIndex = 1;

      for i7=1:length(dim7_idx)
         for i6=1:length(dim6_idx)
            for i5=1:length(dim5_idx)
               for t=1:length(img_idx)

                  #  Position is seeked in bytes. To convert dimension size
                  #  to byte storage size, hdr.dime.bitpix/8 will be
                  #  applied.
                  #
                  # pos = sub2ind([d1 d2 d3 d4 d5 d6 d7], 1, 1, 1, 
			         #             img_idx[t], dim5_idx[i5],dim6_idx[i6],dim7_idx[i7]) -1

                  pos = jl_sub2ind([d1 d2 d3 d4 d5 d6 d7], 1, 1, 1, 
			                      img_idx[t], dim5_idx[i5],dim6_idx[i6],dim7_idx[i7]) - 1
                  
                  #dummy pos
                  pos = 0
                  pos = pos * hdr.dime.bitpix/8

                  if filetype == 2
                     pos2 = convert(Int64, (pos + hdr.dime.vox_offset))
                     seek(io, pos2)
                  else
                     pos2 = convert(Int64, pos)
                     seek(io, pos2)
                  end

                  #  For each frame, fread will read precision of value
                  #  in img_siz times
                  #
                  img[:,currentIndex] = read(io, img_siz) #, sprintf("*%s",precision))
                  currentIndex = currentIndex +1

               end
            end
         end
      end
   end

   #  For complex float32 or complex float64, voxel values
   #  include [real, imag]
   #
   if hdr.dime.datatype == 32 | hdr.dime.datatype == 1792
      img = reshape(img, [2, length(img)/2]);
      img = complex(img(1,:)', img(2,:)');
   end

   end  #fclose(fid);

   #  Update the global min and max values 
   #
   glmax_value = maximum(img[:])
   glmin_value = minimum(img[:])
   # <TODO> How to handle the change of type
   # hdr.dime.glmax = convert(Float64, glmax_value)
   # hdr.dime.glmin = convert(Float64, glmin_value)


   #  old_RGB treat RGB slice by slice, now it is treated voxel by voxel
   #
   old_RGB = true
   if old_RGB == 0
      old_RGB = false
   end
   
   if old_RGB && hdr.dime.datatype == 128 && hdr.dime.bitpix == 24
      # remove squeeze
      img = (reshape(img, [hdr.dime.dim(2:3) 3 hdr.dime.dim(4) length(img_idx) length(dim5_idx) length(dim6_idx) length(dim7_idx)]))
      img = permute(img, [1 2 4 3 5 6 7 8])
   elseif hdr.dime.datatype == 128 & hdr.dime.bitpix == 24
      # remove squeez
      img = (reshape(img, [3 hdr.dime.dim(2:4) length(img_idx) length(dim5_idx) length(dim6_idx) length(dim7_idx)]));
      img = permute(img, [2 3 4 1 5 6 7 8])
   elseif hdr.dime.datatype == 511 & hdr.dime.bitpix == 96
      img = double(img(:));
      img = single((img - min(img))/(max(img) - min(img)))
      # remove squeeze
      indicies = append!([3], hdr.dime.dim(2:4), [length(img_idx)], [length(dim5_idx)], [length(dim6_idx)], [length(dim7_idx)])
      img = (reshape(img, indicies))
      img = permute(img, [2 3 4 1 5 6 7 8])
   else
      # remove squeeze
      indices = append!(hdr.dime.dim[2:4], [length(img_idx)], [length(dim5_idx)], [length(dim6_idx)], [length(dim7_idx)])

      #<TODO> Remove trailing 1's in indices, for now hard coding
      # <TODO> Find a better way later
      ind_tuple = (indices[1], indices[2], indices[3]) #, indices[4], indices[5], indices[6], indices[7])
      img = reshape(img, ind_tuple)
   end

   if ~isempty(img_idx)
      hdr.dime.dim[5] = length(img_idx)
   end

   if ~isempty(dim5_idx)
      hdr.dime.dim[6] = length(dim5_idx)
   end

   if ~isempty(dim6_idx)
      hdr.dime.dim[7] = length(dim6_idx)
   end

   if ~isempty(dim7_idx)
      hdr.dime.dim[8] = length(dim7_idx)
   end

   return [img, hdr]						# read_image
end