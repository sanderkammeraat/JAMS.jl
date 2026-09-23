
@kwdef struct repulsive_soft_disk{T1} <: PairForce
    ontypes::Union{Int64,Vector{Int64}}
    karray::T1
end
function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, force::repulsive_soft_disk)


    if p_i.type in force.ontypes && p_j.type in force.ontypes
        d2R = p_i.R+p_j.R
        if dxn < d2R

            p_i.f += force.karray[get_param_ind(force.ontypes,p_i.type),get_param_ind(force.ontypes,p_j.type)] .* (dxn-d2R) .* dx/dxn
        end
    end
    return p_i

end


@kwdef struct morse{T1, T2}<:PairForce
    ontypes::Union{Int64,Vector{Int64}}
    Dearray::T1
    aarray::T2
    
end


function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, force::morse)

    if p_i.type in force.ontypes && p_j in force.ontypes
        re = p_i.R+p_j.R

        a = force.aarray[get_param_ind(force.ontypes,p_i.type[1]),get_param_ind(force.ontypes,p_j.type[1])]
        De = force.Dearray[get_param_ind(force.ontypes,p_i.type[1]),get_param_ind(force.ontypes,p_j.type[1])]

        p_i.f += -2 * De*a*( exp(-2a*(dxn-re)) - exp(-a*(dxn-re)) ) * dx/dxn

    end
    return p_i

end









function get_param_ind(force_types, particle_type)
    return findfirst(isequal(particle_type),force_types)
end
