

@kwdef struct self_propulsion <: ExternalForce
    ontypes::Union{Int64,Vector{Int64}}
    v0::Float64
end

function contribute_external_force!(p_i,t, dt,rngs_particles, system, force::self_propulsion)
    if p_i.type in force.ontypes
        p_i.f += p_i.zeta .* force.v0 .* p_i.p
    end
    return p_i
end

@kwdef struct planar_rotational_noise <: ExternalForce
    ontypes::Union{Int64,Vector{Int64}}
    Dr::Float64
    normal::SVector{3, Float64} = @SVector [0.0, 0.0, 1.0]
end

function contribute_external_force!(p_i, t, dt, rngs_particles, system, force::planar_rotational_noise)
    if p_i.type in force.ontypes

        η =sqrt( 2*force.Dr ) * rand(rngs_particles[p_i.id],Normal(0, 1))

        p_i.T +=  η .* force.normal .* sqrt(dt)/dt 
    end
    return p_i
end


struct self_align_with_v <: ExternalForce
    ontypes::Union{Int64,Vector{Int64}}
    J::Float64
    unit::Bool
end

function contribute_external_force!(p_i, t, dt,rngs_particles, system, force::self_align_with_v)
    if p_i.type in force.ontypes
        
        if !force.unit
            p_i.T+= force.J*cross(p_i.p,  p_i.v)
        else
            vnorm = norm(p_i.v)
            if vnorm!=0
                p_i.T += force.J*cross(p_i.p,  p_i.v)./vnorm
            end
        end
    end 
    return p_i
end