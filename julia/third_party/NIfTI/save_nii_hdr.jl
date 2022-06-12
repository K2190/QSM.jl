
using Printf

include("data_hist_extra.jl")

#  internal function
  
#  - Jimmy Shen (jimmy@rotman-baycrest.on.ca)

function save_nii_hdr(hdr, io)
   
   # if ~exist("hdr","var") | ~exist("fid","var")
   #    error("Usage: save_nii_hdr(hdr, fid)")
   # end
   
   if hdr.hk.sizeof_hdr != 348
      println("hdr.hk.sizeof_hdr must be 348.")
   end
   
   if hdr.hist.qform_code == 0 && hdr.hist.sform_code == 0
      hdr.hist.sform_code = 1
      hdr.hist.srow_x[1] = hdr.dime.pixdim[2]
      hdr.hist.srow_x[2] = 0
      hdr.hist.srow_x[3] = 0
      hdr.hist.srow_y[1] = 0
      hdr.hist.srow_y[2] = hdr.dime.pixdim[3]
      hdr.hist.srow_y[3] = 0
      hdr.hist.srow_z[1] = 0
      hdr.hist.srow_z[2] = 0
      hdr.hist.srow_z[3] = hdr.dime.pixdim[4]

      hdr_hist_originator = get_dhist_originator();
      hdr.hist.srow_x[4] = (1-hdr_hist_originator[1])*hdr.dime.pixdim[2]
      hdr.hist.srow_y[4] = (1-hdr_hist_originator[2])*hdr.dime.pixdim[3]
      hdr.hist.srow_z[4] = (1-hdr_hist_originator[3])*hdr.dime.pixdim[4]
   end

   write_header(hdr, io)

   return					# save_nii_hdr

end
# %---------------------------------------------------------------------
function write_header(hdr, io)

        #  Original header structures
	#  struct dsr				/* dsr = hdr */
	#       { 
	#       struct header_key hk;            /*   0 +  40       */
	#       struct image_dimension dime;     /*  40 + 108       */
	#       struct data_history hist;        /* 148 + 200       */
	#       };                               /* total= 348 bytes*/
   
   header_key_savenii(io, hdr.hk)
   image_dimension_savenii(io, hdr.dime)
   data_history_savenii(io, hdr.hist)
   
   #  check the file size is 348 bytes
   #
   # fbytes = ftell(fid)
   fbytes = position(io)
   if fbytes != 348
      msg = @sprintf("Header size is not 348 bytes.")
      println(msg)
   end
    
   return					# write_header

end

# %---------------------------------------------------------------------
function header_key_savenii(io, hk)
   
   seek(io,0)

	#  Original header structures    
	#  struct header_key                      /* header key      */ 
	#       {                                /* off + size      */
	#       int sizeof_hdr                   /*  0 +  4         */
	#       char data_type[10];              /*  4 + 10         */
	#       char db_name[18];                /* 14 + 18         */
	#       int extents;                     /* 32 +  4         */
	#       short int session_error;         /* 36 +  2         */
	#       char regular;                    /* 38 +  1         */
	#       char dim_info;   # char hkey_un0;        /* 39 +  1 */
	#       };                               /* total=40 bytes  */
        
   write(io, hk.sizeof_hdr[1])	# must be 348.

   # data_type = sprintf("%-10s",hk.data_type);	# ensure it is 10 chars from left
   # fwrite(fid, data_type(1:10), "uchar");
   pad = zeros(1, 10-length(hk.data_type))

   #<TODO> Check this in debugging
   #hk.data_type = hk.data_type  #  char(pad)]
   # fwrite(fid, hk.data_type[1:10], "uchar")

   hk_length = length(hk.data_type)
   for i=hk_length+1:10
      append!(hk.data_type, 0)
   end
   write(io, hk.data_type[1:10])
    
   # db_name   = sprintf("%-18s", hk.db_name);	# ensure it is 18 chars from left
   # fwrite(fid, db_name(1:18), "uchar");
   # pad = zeros(1, 18-length(hk.db_name))
   # hk.db_name = [hk.db_name  char(pad)]
   # fwrite(fid, hk.db_name[1:18], "uchar")
   hk_length = length(hk.db_name)
   for i=hk_length+1:18
      append!(hk.db_name, 0)
   end
   write(io, hk.db_name[1:18])

   # fwrite(fid, hk.extents[1],       "int32")
   # fwrite(fid, hk.session_error[1], "int16")
   # fwrite(fid, hk.regular[1],       "uchar")	          # might be uint8
   write(io, hk.extents[1])
   write(io, hk.session_error[1])
   write(io, hk.regular[1])	          # might be uint8   
    
   # fwrite(fid, hk.hkey_un0[1],    "uchar")
   # fwrite(fid, hk.hkey_un0[1],    "uint8")

   # fwrite(fid, hk.dim_info[1],      "uchar")
   write(io, hk.dim_info[1])
   
   return					# header_key
end

# %---------------------------------------------------------------------
function image_dimension_savenii(io, dime)

	#  Original header structures        
	#  struct image_dimension
	#       {                                /* off + size      */
	#       short int dim[8];                /* 0 + 16          */
	#       float intent_p1;   # char vox_units[4];   /* 16 + 4       */
	#       float intent_p2;   # char cal_units[8];   /* 20 + 4       */
	#       float intent_p3;   % char cal_units[8];   /* 24 + 4       */
	#       short int intent_code;   # short int unused1;   /* 28 + 2 */
	#       short int datatype;              /* 30 + 2          */
	#       short int bitpix;                /* 32 + 2          */
	#       short int slice_start;   # short int dim_un0;   /* 34 + 2 */
	#       float pixdim[8];                 /* 36 + 32         */
	#			/*
	#				pixdim[] specifies the voxel dimensions:
	#				pixdim[1] - voxel width
	#				pixdim[2] - voxel height
	#				pixdim[3] - interslice distance
	#				pixdim[4] - volume timing, in msec
	#					..etc
	#			*/
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
	
   write(io, dime.dim[1:8])
   write(io, dime.intent_p1[1])   # ,  "float32")
   write(io, dime.intent_p2[1])   #,  "float32")
   write(io, dime.intent_p3[1])   #,  "float32")
   write(io, dime.intent_code[1])   #,  "int16")
   write(io, dime.datatype[1])   #,     "int16")
   write(io, dime.bitpix[1])   #,       "int16")
   write(io, dime.slice_start[1])   #,  "int16")
   write(io, dime.pixdim[1:8])   #,   "float32")
   write(io, dime.vox_offset[1])   #, "float32")
   write(io, dime.scl_slope[1])   #,  "float32")
   write(io, dime.scl_inter[1])   #,  "float32")
   write(io, dime.slice_end[1])   #,    "int16")
   write(io, dime.slice_code[1])   #,   "uchar")
   write(io, dime.xyzt_units[1])   #,   "uchar")
   write(io, dime.cal_max[1])   #,    "float32")
   write(io, dime.cal_min[1])   #,    "float32")
   write(io, dime.slice_duration[1])   #, "float32")
   write(io, dime.toffset[1])   #,    "float32")
   write(io, dime.glmax[1])       #,        "int32")
   write(io, dime.glmin[1])   #,        "int32")
   
   return;					# image_dimension
end

# %---------------------------------------------------------------------
function data_history_savenii(io, hist)
    
	# Original header structures
	#struct data_history       
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
	
   # descrip     = sprintf("%-80s", hist.descrip)     # 80 chars from left
   # fwrite(fid, descrip(1:80),    "uchar")

   # pad = zeros(1, 80-length(hist.descrip))
   # hist.descrip = [hist.descrip  char(pad)]
   descrip_length = length(hist.descrip)
   for i=descrip_length+1:80
      append!(hist.descrip, 0)
   end

   write(io, hist.descrip[1:80])   #, "uchar")
    
   # aux_file    = sprintf("%-24s", hist.aux_file);    # 24 chars from left
   # fwrite(fid, aux_file(1:24),   "uchar")
   # pad = zeros(1, 24-length(hist.aux_file))
   # hist.aux_file = [hist.aux_file  char(pad)]
   aux_file_length = length(hist.aux_file)
   for i=aux_file_length+1:24
      append!(hist.aux_file, 0)
   end
   write(io, hist.aux_file[1:24])   #, "uchar")
    
   write(io, hist.qform_code)   #,,    "int16")
   write(io, hist.sform_code)   #,,    "int16")
   write(io, hist.quatern_b)   #,,   "float32")
   write(io, hist.quatern_c)   #,,   "float32")
   write(io, hist.quatern_d)   #,,   "float32")
   write(io, hist.qoffset_x)   #,,   "float32")
   write(io, hist.qoffset_y)   #,,   "float32")
   write(io, hist.qoffset_z)   #,,   "float32")
   write(io, hist.srow_x[1:4])   #,, "float32")
   write(io, hist.srow_y[1:4])   #,, "float32")
   write(io, hist.srow_z[1:4])   #,, "float32")

   # intent_name = sprintf("%-16s", hist.intent_name)	# 16 chars from left
   # fwrite(fid, intent_name(1:16),    "uchar")
   # pad = zeros(1, 16-length(hist.intent_name))
   # hist.intent_name = [hist.intent_name  char(pad)]
   hist_length = length(hist.intent_name)
   for i=hist_length+1:16
      append!(hist.intent_name, 0)
   end
   write(io, hist.intent_name[1:16])     #, "uchar")
    
   # magic	= sprintf("%-4s", hist.magic)		# 4 chars from left
   # fwrite(fid, magic(1:4),           "uchar")
   # pad = zeros(1, 4-length(hist.magic))
   # hist.magic = [hist.magic  char(pad)]
   hist_length = length(hist.magic)
   for i=hist_length+1:4
      append!(hist.magic, 0)
   end
   write(io, hist.magic[1:4])    #,        "uchar")
    
   return					# data_history

end