module Save


function polar_particle!(current_frame_group, current_particle_state, current_field_state, n, Tsave, t,framecounter)

    current_frame_group["n"] = n

    current_frame_group["t"] = t

    current_frame_group["id"] = current_particle_state.id

    current_frame_group["type"] = current_particle_state.type

    current_frame_group["R"] = current_particle_state.R

    current_frame_group["x"]   = reinterpret(reshape, Float64, current_particle_state.x)
    current_frame_group["xuw"] = reinterpret(reshape, Float64, current_particle_state.xuw)
    current_frame_group["v"]   = reinterpret(reshape, Float64, current_particle_state.v)
    current_frame_group["p"]   = reinterpret(reshape, Float64, current_particle_state.p)
    current_frame_group["q"]   = reinterpret(reshape, Float64, current_particle_state.q)

    return current_frame_group


end


end