mutable struct data_history_extra
    originator::Vector{Int16}
    flip_orient::Vector{UInt8}
    rot_orient::Vector{UInt8}
  end

global dhist_extra = data_history_extra(zeros(3), zeros(3), zeros(3))

function get_dhist_originator()
    return dhist_extra.originator
end

function get_dhist_flip_orient()
    return dhist_extra.originator
end

function get_dhist_rot_orient()
    return dhist_extra.rot_orient
end

function set_dhist_originator(originator)
    global dhist_extra.originator[1] = originator[1]
    global dhist_extra.originator[2] = originator[2]
    global dhist_extra.originator[3] = originator[3]
end

function set_dhist_flip_orient(flip_orient)
    global dhist_extra.flip_orient[1] = flip_orient[1]
    global dhist_extra.flip_orient[2] = flip_orient[2]
    global dhist_extra.flip_orient[3] = flip_orient[3]
end

function set_dhist_rot_orient(rot_orient)
    global dhist_extra.rot_orient[1] = rot_orient[1]
    global dhist_extra.rot_orient[2] = rot_orient[2]
    global dhist_extra.rot_orient[3] = rot_orient[3]
end

