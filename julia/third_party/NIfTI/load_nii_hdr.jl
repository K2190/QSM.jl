using Printf

mutable struct data_history
	descrip::Vector{UInt8}         #  /* 0 + 80          */
	aux_file::Vector{UInt8}        #  /* 80 + 24         */
	qform_code::Int16              #  /* 104 + 2         */
	sform_code::Int16              #  /* 106 + 2         */
	quatern_b::Float32             #  /* 108 + 4         */
	quatern_c::Float32             #  /* 112 + 4         */
	quatern_d::Float32             #  /* 116 + 4         */
	qoffset_x::Float32             #  /* 120 + 4         */
	qoffset_y::Float32             #  /* 124 + 4         */
	qoffset_z::Float32             #  /* 128 + 4         */
	srow_x::Vector{Float32}        #  /* 132 + 16        */
	srow_y::Vector{Float32}        #  /* 148 + 16        */
	srow_z::Vector{Float32}        #  /* 164 + 16        */
	intent_name::Vector{UInt8}     #  /* 180 + 16        */
	magic::Vector{UInt8}           #  /* 196 + 4         */
 end                             #  /* total=200 bytes */
    
 mutable struct image_dimension
   dim::Vector{Int16}                #/* 0 + 16          */
   intent_p1::Float32                # char vox_units[4];   /* 16 + 4       */ 
   intent_p2::Float32                # char cal_units[8];   /* 20 + 4       */
   intent_p3::Float32                # char cal_units[8];   /* 24 + 4       */
   intent_code::Int16                # short int unused1;   /* 28 + 2 */
   datatype::Int16                   # /* 30 + 2          */
   bitpix::Int16                     # /* 32 + 2          */
   slice_start::Int16                # short int dim_un0;   /* 34 + 2 */
   pixdim::Vector{Float32}                #  /* 36 + 32         */
   vox_offset::Float32               # /* 68 + 4          */
   scl_slope::Float32                # float roi_scale;     /* 72 + 4 */
   scl_inter::Float32                # float funused1;      /* 76 + 4 */
   slice_end::Float16                # float funused2;      /* 80 + 2 */
   slice_code::Int8                  # float funused2;      /* 82 + 1 */
   xyzt_units::Int8                  # float funused2;      /* 83 + 1 */
   cal_max::Float32                  #       /* 84 + 4          */
   cal_min::Float32                  #       /* 88 + 4          */
   slice_duration::Float32           # int compressed; /* 92 + 4 */
   toffset::Float32                  # int verified;          /* 96 + 4 */
   glmax::Int32                      #/* 100 + 4         */
   glmin::Int32                      # /* 104 + 4         */
end                                  #/* total=108 bytes */


mutable struct header_key          # /* header key      */ 
   sizeof_hdr::Int32               # /*  0 +  4         */
   data_type::Vector{Int8}         # /*  4 + 10         */
   db_name::Vector{Int8}           # /* 14 + 18         */
   extents::Int32                  # /* 32 +  4         */
   session_error::Int16            # /* 36 +  2         */
   regular::Int8                   # /* 38 +  1         */
   dim_info::Int8                  # /* 39 +  1 */
end  

mutable struct dsr
  hk::header_key             # /*   0 +  40       */
  dime::image_dimension      #/ *  40 + 108       */
  hist::data_history         # /* 148 + 200       */
end 

#  internal function
#
#  - Jimmy Shen (jimmy@rotman-baycrest.on.ca)

function load_nii_hdr(fileprefix)

   # if ~exist("fileprefix","var"),
   #    error("Usage: [hdr, filetype, fileprefix, machine] = load_nii_hdr(filename)");
   # end

   
   machine = "ieee-le"
   new_ext = 0

   if occursin(".nii",fileprefix) && fileprefix[end-3:end] == ".nii"
      new_ext = 1
      #fileprefix[end-3:end] = ""
      fileprefix2 = fileprefix[1:end-4] 
      fileprefix = fileprefix2
   end

   if occursin(".hdr",fileprefix) && fileprefix[end-3:end] == ".hdr"
      #fileprefix[end-3:end] = ""
      fileprefix2 = fileprefix[1:end-4] 
      fileprefix = fileprefix2
   end

   if occursin(".img",fileprefix) && fileprefix[end-3:end] == ".img"
      # fileprefix[end-3:end] = ""
      fileprefix2 = fileprefix[1:end-4] 
      fileprefix = fileprefix2
   end

   if new_ext > 0
      fn = @sprintf("%s.nii",fileprefix)

      if isfile(fn) != 1
         msg = @sprintf("Cannot find file %s.nii ", fileprefix)
         println(msg)
      end
   else
      fn = @sprintf("%s.hdr",fileprefix)

      if isfile(fn) != 1
         msg = @sprintf("Cannot find file %s.hdr", fileprefix)
         println(msg)
      end
   end

   #fid = open(fn,"r")  #<TODO> ,machine)
   
   fn_data = zeros(UInt32, 1)
   hdr_read = 0
   open(fn) do io
      #  first try reading the little endian 
      read!(io, fn_data)

      if fn_data[1] == 348
         global hdr = read_header(io, 0)
      else
         #  Now reading the big endian 
         fn_data .= bswap.(fn_data)

         if fn_data[1] == 348
            
            global hdr = read_header(io, 1)
         end
      end
   end

  
   hist_magic = hdr.hist.magic
   indecies = findall(x->x==0, hist_magic)
   indx = indecies[1]
   hist_magic = hist_magic[1 : indx-1]
   hist_magic = convert(Vector{UInt8}, hist_magic)
   hist_magic_str = String(hist_magic)

   if hist_magic_str == "n+1"
      filetype = 2
   elseif hist_magic_str == "ni1"
      filetype = 1
   else
      filetype = 0
   end

   return [hdr, filetype, fileprefix] 

end

#---------------------------------------------------------------------
function read_header(io, swap)

        #  Original header structures
	#  struct dsr
	#       { 
	#       struct header_key hk;            /*   0 +  40       */
	#       struct image_dimension dime;     /*  40 + 108       */
	#       struct data_history hist;        /* 148 + 200       */
	#       };                               /* total= 348 bytes*/

   
    dsr_hk   = header_key(io, swap)
    dsr_dime = image_dimension(io, swap)
    dsr_hist = data_history(io, swap)

    dsr_obj = dsr(dsr_hk, dsr_dime, dsr_hist)

   #  #  For Analyze data format
   #  #
   #  if dsr_obj.hist.magic != "n+1" && dsr_obj.hist.magic != "ni1"
   #    dsr_obj.hist.qform_code = 0
   #    dsr_obj.hist.sform_code = 0
   #  end

   #  return dsr_obj					# read_header
end

#---------------------------------------------------------------------
function header_key(io, swap)

   # fseek(fid,0,"bof");
   seek(io, 0)

	#  Original header structures	
	#  struct header_key                     /* header key      */ 
	#       {                                /* off + size      */
	#       int sizeof_hdr                   /*  0 +  4         */
	#       char data_type[10];              /*  4 + 10         */
	#       char db_name[18];                /* 14 + 18         */
	#       int extents;                     /* 32 +  4         */
	#       short int session_error;         /* 36 +  2         */
	#       char regular;                    /* 38 +  1         */
	#       char dim_info;   # char hkey_un0;        /* 39 +  1 */
	#       };                               /* total=40 bytes  */
	#
	# int sizeof_header   Should be 348.
	# char regular        Must be "r" to indicate that all images and 
	#                     volumes are the same size. 


    # NOT Needed in JULIA
    #v6 = version;
    #if str2num(v6(1))<6
    #   directchar = "*char"
    #else
    #   directchar = "uchar=>char"
    #end

  

      int32_data = zeros(Int32, 1)
      int16_data = zeros(Int16, 1)
      char_data_10 = zeros(UInt8, 10)
      char_data_18 = zeros(UInt8, 18)
      char_data = zeros(UInt8, 1)

      read!(io, int32_data)
      if swap == 1
        int32_data .= bswap.(int32_data)
      end
      hk_sizeof_hdr = int32_data[1]
      #Should be 348
      if hk_sizeof_hdr != 348
         println("ERROR: Data incorrect")
      end

      read!(io, char_data_10)
      if swap == 1
        char_data_10 .= bswap.(char_data_10)
      end
      hk_data_type = char_data_10

      read!(io, char_data_18)
      if swap == 1
        char_data_18 .= bswap.(char_data_18)
      end
      hk_db_name = char_data_18

      read!(io, int32_data)
      if swap == 1
        int32_data .= bswap.(int32_data)
      end
      hk_extents = int32_data[1]

      read!(io, int16_data)
      if swap == 1
        int16_data .= bswap.(int16_data)
      end
      hk_session_error = int16_data[1]

      read!(io, char_data)
      if swap == 1
        char_data .= bswap.(char_data)
      end
      hk_regular = char_data[1]

      read!(io, char_data)
      if swap == 1
        char_data .= bswap.(char_data)
      end
      hk_dim_info = char_data[1]

      h_key = header_key(hk_sizeof_hdr, hk_data_type, hk_db_name, hk_extents, 
                          hk_session_error, hk_regular, hk_dim_info)
      return	h_key 


   end

#---------------------------------------------------------------------
function image_dimension(io, swap)

	#  Original header structures    
	#  struct image_dimension
	#       {                                /* off + size      */
	#       short int dim[8];                /* 0 + 16          */
        #       /*
        #           dim[0]      Number of dimensions in database; usually 4. 
        #           dim[1]      Image X dimension;  number of *pixels* in an image row. 
        #           dim[2]      Image Y dimension;  number of *pixel rows* in slice. 
        #           dim[3]      Volume Z dimension; number of *slices* in a volume. 
        #           dim[4]      Time points; number of volumes in database
        #       */
	#       float intent_p1;   # char vox_units[4];   /* 16 + 4       */
	#       float intent_p2;   # char cal_units[8];   /* 20 + 4       */
	#       float intent_p3;   # char cal_units[8];   /* 24 + 4       */
	#       short int intent_code;   # short int unused1;   /* 28 + 2 */
	#       short int datatype;              /* 30 + 2          */
	#       short int bitpix;                /* 32 + 2          */
	#       short int slice_start;   # short int dim_un0;   /* 34 + 2 */
	#       float pixdim[8];                 /* 36 + 32         */
	#	/*
	#		pixdim[] specifies the voxel dimensions:
	#		pixdim[1] - voxel width, mm
	#		pixdim[2] - voxel height, mm
	#		pixdim[3] - slice thickness, mm
	#		pixdim[4] - volume timing, in msec
	#					..etc
	#	*/
	#       float vox_offset;                /* 68 + 4          */
	#       float scl_slope;   # float roi_scale;     /* 72 + 4 */
	#       float scl_inter;   # float funused1;      /* 76 + 4 */
	#       short slice_end;   # float funused2;      /* 80 + 2 */
	#       char slice_code;   # float funused2;      /* 82 + 1 */
	#       char xyzt_units;   # float funused2;      /* 83 + 1 */
	#       float cal_max;                   /* 84 + 4          */
	#       float cal_min;                   /* 88 + 4          */
	#       float slice_duration;   # int compressed; /* 92 + 4 */
	#       float toffset;   # int verified;          /* 96 + 4 */
	#       int glmax;                       /* 100 + 4         */
	#       int glmin;                       /* 104 + 4         */
	#       };                               /* total=108 bytes */
	
   
   int16_data = zeros(Int16, 1)
   int16_data_8 = zeros(Int16, 8)
   int32_data = zeros(Int32, 1)
   uchar_data = zeros(Int8, 1)   
   float32_data = zeros(Float32, 1)
   float32_data_8 = zeros(Float32, 8)


    #dime.dim        = fread(fid,8,"int16")'
    read!(io, int16_data_8)
    if swap == 1
      int16_data_8 .= bswap.(int16_data_8)
    end
    dime_dim = int16_data_8

   #  println("IN IMG DIMENSION")
   #  println(dime_dim)

    #dime.intent_p1  = fread(fid,1,"float32")'
    read!(io, float32_data)
    if swap == 1
      float32_data .= bswap.(float32_data)
    end
    dime_intent_p1 = float32_data[1]

    #dime.intent_p2  = fread(fid,1,"float32")'
    read!(io, float32_data)
    if swap == 1
      float32_data .= bswap.(float32_data)
    end
    dime_intent_p2 = float32_data[1]
    
    #dime.intent_p3  = fread(fid,1,"float32")'
    read!(io, float32_data)
    if swap == 1
      float32_data .= bswap.(float32_data)
    end
    dime_intent_p3 = float32_data[1]
    
    #dime.intent_code = fread(fid,1,"int16")'
    read!(io, int16_data)
    if swap == 1
      int16_data .= bswap.(int16_data)
    end
    dime_intent_code = int16_data[1]
    
    #dime.datatype   = fread(fid,1,"int16")'
    read!(io, int16_data)
    if swap == 1
      int16_data .= bswap.(int16_data)
    end
    dime_datatype = int16_data[1]
    
    #dime.bitpix     = fread(fid,1,"int16")'
    read!(io, int16_data)
    if swap == 1
      int16_data .= bswap.(int16_data)
    end
    dime_bitpix = int16_data[1]
    
    #dime.slice_start = fread(fid,1,"int16")'
    read!(io, int16_data)
    if swap == 1
      int16_data .= bswap.(int16_data)
    end
    dime_slice_start = int16_data[1]
    
    #dime.pixdim     = fread(fid,8,"float32")'
    read!(io, float32_data_8)
    if swap == 1
      float32_data_8 .= bswap.(float32_data_8)
    end
    dime_pixdim = float32_data_8
    
    #dime.vox_offset = fread(fid,1,"float32")'
    read!(io, float32_data)
    if swap == 1
      float32_data .= bswap.(float32_data)
    end
    dime_vox_offset = float32_data[1]
    
    #dime.scl_slope  = fread(fid,1,"float32")'
    read!(io, float32_data)
    if swap == 1
      float32_data .= bswap.(float32_data)
    end
    dime_scl_slope = float32_data[1]
    
    #dime.scl_inter  = fread(fid,1,"float32")'
    read!(io, float32_data)
    if swap == 1
      float32_data .= bswap.(float32_data)
    end
    dime_scl_inter = float32_data[1]
    
    #dime.slice_end  = fread(fid,1,"int16")'
    read!(io, int16_data)
    if swap == 1
      int16_data .= bswap.(int16_data)
    end
    dime_slice_end = int16_data[1]
    
    #dime.slice_code = fread(fid,1,"uchar")'
    read!(io, uchar_data)
    if swap == 1
      int16_data .= bswap.(uchar_data)
    end
    dime_slice_code = uchar_data[1]
    
    #dime.xyzt_units = fread(fid,1,"uchar")'
    read!(io, uchar_data)
    if swap == 1
      uchar_data .= bswap.(uchar_data)
    end
    dime_xyzt_units = uchar_data[1]
    
    #dime.cal_max    = fread(fid,1,"float32")'
    read!(io, float32_data)
    if swap == 1
      float32_data .= bswap.(float32_data)
    end
    dime_cal_max = float32_data[1]
    
    #dime.cal_min    = fread(fid,1,"float32")'
    read!(io, float32_data)
    if swap == 1
      float32_data .= bswap.(float32_data)
    end
    dime_cal_min = float32_data[1]
    
    #dime.slice_duration = fread(fid,1,"float32")'
    read!(io, float32_data)
    if swap == 1
      float32_data .= bswap.(float32_data)
    end
    dime_slice_duration = float32_data[1]
    
    #dime.toffset    = fread(fid,1,"float32")'
    read!(io, float32_data)
    if swap == 1
      float32_data .= bswap.(float32_data)
    end
    dime_toffset = float32_data[1]
    
    #dime.glmax      = fread(fid,1,"int32")'
    read!(io, int32_data)
    if swap == 1
      int32_data .= bswap.(int32_data)
    end
    dime_glmax = int32_data[1]
    
    #dime.glmin      = fread(fid,1,"int32")'
    read!(io, int32_data)
    if swap == 1
      int32_data .= bswap.(int32_data)
    end
    dime_glmin = int32_data[1]
    
    dime = image_dimension(dime_dim, dime_intent_p1, dime_intent_p2, dime_intent_p3, dime_intent_code,
                           dime_datatype, dime_bitpix, dime_slice_start, dime_pixdim, dime_vox_offset,
                           dime_scl_slope, dime_scl_inter, dime_slice_end, dime_slice_code, 
                           dime_xyzt_units, dime_cal_max, dime_cal_min, dime_slice_duration,
                           dime_toffset, dime_glmax,dime_glmin)

   
        
    return dime 				# image_dimension



end

#---------------------------------------------------------------------
function data_history(io, swap)
        
	#  Original header structures
	#  struct data_history       
	#       {                                /* off + size      */
	#       char descrip[80];                /* 0 + 80          */
	#       char aux_file[24];               /* 80 + 24         */
	#       short int qform_code;            /* 104 + 2         */
	#       short int sform_code;            /* 106 + 2         */
	#       float quatern_b;                 /* 108 + 4         */
	#       float quatern_c;                 /* 112 + 4         */
	#       float quatern_d;                 /* 116 + 4         */
	#       float qoffset_x;                 /* 120 + 4         */
	#       float qoffset_y;                 /* 124 + 4         */
	#       float qoffset_z;                 /* 128 + 4         */
	#       float srow_x[4];                 /* 132 + 16        */
	#       float srow_y[4];                 /* 148 + 16        */
	#       float srow_z[4];                 /* 164 + 16        */
	#       char intent_name[16];            /* 180 + 16        */
	#       char magic[4];   # int smin;     /* 196 + 4         */
	#       };                               /* total=200 bytes */

   #  v6 = version;
   #  if str2num(v6(1))<6
   #     directchar = "*char";
   #  else
   #     directchar = "uchar=>char";
   #  end

    int16_data = zeros(Int16, 1)
    int16_data_orgn = zeros(Int16, 5)
    uchar_data80 = zeros(UInt8, 80) 
    uchar_data24 = zeros(UInt8, 24) 
    uchar_data16 = zeros(UInt8, 16) 
    uchar_data4 = zeros(UInt8, 4) 
    float32_data = zeros(Float32, 1)
    float32_data_4_x = zeros(Float32, 4)
    float32_data_4_y = zeros(Float32, 4)
    float32_data_4_z = zeros(Float32, 4)

        #hist.descrip     = deblank(fread(fid,80,directchar)')
        read!(io, uchar_data80)
        if swap == 1
          uchar_data80 .= bswap.(uchar_data80)
        end
        hist_descrip = uchar_data80
        
        #hist.aux_file    = deblank(fread(fid,24,directchar)')
        read!(io, uchar_data24)
        if swap == 1
          uchar_data24 .= bswap.(uchar_data24)
        end
        hist_aux_file = uchar_data24
        
        #hist.qform_code  = fread(fid,1,"int16")'
        read!(io, int16_data)
        if swap == 1
          int16_data .= bswap.(int16_data)
        end
        hist_qform_code = int16_data[1]
        
        #hist.sform_code  = fread(fid,1,"int16")'
        read!(io, int16_data)
        if swap == 1
          int16_data .= bswap.(int16_data)
        end
        hist_sform_code = int16_data[1]
        
        #hist.quatern_b   = fread(fid,1,"float32")'
        read!(io, float32_data)
        if swap == 1
          float32_data .= bswap.(float32_data)
        end
        hist_quatern_b = float32_data[1]
        
        #hist.quatern_c   = fread(fid,1,"float32")'
        read!(io, float32_data)
        if swap == 1
          float32_data .= bswap.(float32_data)
        end
        hist_quatern_c = float32_data[1]
        
        #hist.quatern_d   = fread(fid,1,"float32")'
        read!(io, float32_data)
        if swap == 1
          float32_data .= bswap.(float32_data)
        end
        hist_quatern_d = float32_data[1]
        
        #hist.qoffset_x   = fread(fid,1,"float32")'
        read!(io, float32_data)
        if swap == 1
          float32_data .= bswap.(float32_data)
        end
        hist_qoffset_x = float32_data[1]
        
        #hist.qoffset_y   = fread(fid,1,"float32")'
        read!(io, float32_data)
        if swap == 1
          float32_data .= bswap.(float32_data)
        end
        hist_qoffset_y = float32_data[1]
        
        #hist.qoffset_z   = fread(fid,1,"float32")'
        read!(io, float32_data)
        if swap == 1
          float32_data .= bswap.(float32_data)
        end
        hist_qoffset_z = float32_data[1]
        
        #hist.srow_x      = fread(fid,4,"float32")'
        read!(io, float32_data_4_x)
        if swap == 1
          float32_data_4_x .= bswap.(float32_data_4_x)
        end
        hist_srow_x = float32_data_4_x
        
        #hist.srow_y      = fread(fid,4,"float32")'
        read!(io, float32_data_4_y)
        if swap == 1
          float32_data_4_y .= bswap.(float32_data_4_y)
        end
        hist_srow_y = float32_data_4_y
        
        #hist.srow_z      = fread(fid,4,"float32")'
        read!(io, float32_data_4_z)
        if swap == 1
          float32_data_4_z .= bswap.(float32_data_4_z)
        end
        hist_srow_z = float32_data_4_z
        
        #hist.intent_name = deblank(fread(fid,16,directchar)')
        read!(io, uchar_data16)
        if swap == 1
          uchar_data16 .= bswap.(uchar_data16)
        end
        hist_intent_name = uchar_data16
        
        #hist.magic       = deblank(fread(fid,4,directchar)')
        read!(io, uchar_data4)
        if swap == 1
          uchar_data4 .= bswap.(uchar_data4)
        end
        hist_magic = uchar_data4

        seek(io,253)
        read!(io, int16_data_orgn)
        if swap == 1
          int16_data_orgn .= bswap.(int16_data_orgn)
        end
        hist_originator = int16_data_orgn

        # Initialze to 0.0, will be filled later
        flip_orient= zeros(UInt8, 3) 
        rot_orient= zeros(UInt8, 3) 
      

      hist = data_history(hist_descrip,  hist_aux_file,  hist_qform_code,  hist_sform_code,  hist_quatern_b,
                            hist_quatern_c,  hist_quatern_d,  hist_qoffset_x,  hist_qoffset_y,  hist_qoffset_z,
                            hist_srow_x,  hist_srow_y,  hist_srow_z,  hist_intent_name,  hist_magic)
                           #  hist_originator, flip_orient, rot_orient)

      return hist 


end
