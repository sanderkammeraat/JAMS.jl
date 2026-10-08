

"""
    self_propulsion(ontypes, v0)
    self_propulsion(; ontypes, v0)

Constant self-propulsion force along the polarity vector:
`f += zeta * v0 * p`. For an overdamped particle this gives a self-propulsion speed `v0`.

# Fields

- `ontypes`: particle type (`Int`) or types (`Vector{Int}`) this force acts on
- `v0`: self-propulsion speed

Requires particle fields `type`, `zeta`, `p` and `f`.
"""
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

"""
    planar_rotational_noise(; ontypes, Dr, normal=[0, 0, 1])

Gaussian white-noise torque about the axis `normal`:
`T += sqrt(2Dr) * η  sqrt(dt) * normal`, with `η` a standard normal random number drawn
from the particle's own random number generator every timestep.

Combined with [`DOFevolvers.overdamped_pqT_xyc`](@ref JAMS.DOFevolvers.overdamped_pqT_xyc), the polarity angle in the xy-plane
then diffuses with rotational diffusion coefficient `Dr / zeta_R^2` (so `Dr` for the default
`zeta_R = 1`).

# Fields

- `ontypes`: particle type (`Int`) or types (`Vector{Int}`) this force acts on
- `Dr`: noise strength (rotational diffusion coefficient for `zeta_R = 1`)
- `normal`: rotation axis (**Default**: `[0, 0, 1]`, rotations in the xy-plane)

Requires particle fields `id`, `type` and `T`.
"""
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

@kwdef struct translational_noise <: ExternalForce
    ontypes::Union{Int64,Vector{Int64}}
    Dtarray::SVector{3, Float64} 
end


#continue here...
function contribute_external_force!(p_i, t, dt, rngs_particles, system, force::translational_noise)
    if p_i.type in force.ontypes

        eta_x = rand(rngs_particles[p_i.id],Normal(0, 1))*sqrt(2*force.Dtarray[1])* sqrt(dt)/dt 
        eta_y = rand(rngs_particles[p_i.id],Normal(0, 1))*sqrt(2*force.Dtarray[2])* sqrt(dt)/dt 
        eta_z = rand(rngs_particles[p_i.id],Normal(0, 1))*sqrt(2*force.Dtarray[3])* sqrt(dt)/dt 
        p_i.f +=  SVector{3, Float64}(eta_x, eta_y, eta_z)
    end
    return p_i
end


"""
    self_align_with_v(ontypes, J, unit)

Torque that rotates the polarity `p` towards the particle's own velocity `v`
(self-alignment): `T += J * (p × v)`, or `T += J * (p × v / |v|)` if `unit` is `true`.
The velocity used is the one set by the DOF evolver in the previous timestep, so this force
does nothing in the first timestep if the particles start with `v = 0`.

# Fields

- `ontypes`: particle type (`Int`) or types (`Vector{Int}`) this force acts on
- `J`: alignment strength
- `unit`: if `true`, use the direction of `v` only (the torque does not grow with speed)

Requires particle fields `type`, `p`, `v` and `T`.
"""
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