
export repulsive_soft_disk
@kwdef struct repulsive_soft_disk{T1} <: PairForce
    ontypes::Union{Int64,Vector{Int64}}
    karray::T1
end
function contribute_pair_force!(i,p_i, p_j , current_particle_state, dx, dxn, t, dt,rngs_particles, system, force::repulsive_soft_disk)


    if p_i.type in force.ontypes && p_j.type in force.ontypes
    d2R = p_i.R+p_j.R
        if dxn < d2R

            current_particle_state.f[i]+= force.karray[get_param_ind(force.ontypes,p_i.type),get_param_ind(force.ontypes,p_j.type)] .* (dxn-d2R) .* dx/dxn
        end
    end

end
function get_param_ind(force_types, particle_type)
    return findfirst(isequal(particle_type),force_types)
end
