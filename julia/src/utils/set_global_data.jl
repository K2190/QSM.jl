

global cur_nii_hdr_dime_datatype = 0

function set_nii_obj_info(dime_datatype)
    global cur_nii_hdr_dime_datatype = dime_datatype
end

function get_nii_obj_dime_datatype()
    return cur_nii_hdr_dime_datatype
end