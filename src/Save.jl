module Save


function polar_particle!(current_frame_group, current_particle_state, current_field_state, n, Tsave, t,framecounter)

    current_frame_group["n"] = n

    current_frame_group["t"] = t

    current_frame_group["id"] = current_particle_state.id

    current_frame_group["type"] = current_particle_state.type

    current_frame_group["R"] = current_particle_state.R

    current_frame_group["x"] = [p_i.x[1] for p_i in current_particle_state]
    current_frame_group["y"] = [p_i.x[2] for p_i in current_particle_state]
    current_frame_group["z"] = [p_i.x[3] for p_i in current_particle_state]

    current_frame_group["xuw"] = [p_i.xuw[1] for p_i in current_particle_state]
    current_frame_group["yuw"] = [p_i.xuw[2] for p_i in current_particle_state]
    current_frame_group["zuw"] = [p_i.xuw[3] for p_i in current_particle_state]

    current_frame_group["vx"] = [p_i.v[1] for p_i in current_particle_state]
    current_frame_group["vy"] = [p_i.v[2] for p_i in current_particle_state]
    current_frame_group["vz"] = [p_i.v[3] for p_i in current_particle_state]

    current_frame_group["px"] = [p_i.p[1] for p_i in current_particle_state]
    current_frame_group["py"] = [p_i.p[2] for p_i in current_particle_state]
    current_frame_group["pz"] = [p_i.p[3] for p_i in current_particle_state]

    current_frame_group["qx"] = [p_i.q[1] for p_i in current_particle_state]
    current_frame_group["qy"] = [p_i.q[2] for p_i in current_particle_state]
    current_frame_group["qz"] = [p_i.q[3] for p_i in current_particle_state]

    return current_frame_group


end


end