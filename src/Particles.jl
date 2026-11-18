module Particles

using StaticArrays
using LinearAlgebra
using Distributions
using Random
using SparseArrays

"""
    Particle

Abstract supertype of all particle types. A concrete particle type is an immutable struct
whose fields are the degrees of freedom and parameters of one particle. The simulation
stores all particles in a [`ParticleState`](@ref JAMS.ParticleState) (a `StructArray`), so every field becomes
a column that forces and DOF evolvers can read and update.

A particle type needs at least the fields `id`, `type`, `x`, `xuw`, `f` and `T`; see the
*Extending JAMS* page for details.
"""
abstract type Particle end

"""
    RigidBody

Abstract supertype for extended (rigid-body) particles.
"""
abstract type RigidBody end

"""
    Polar(; id, type, x, p, kwargs...)

Spherical particle with a position and a polarity (orientation) vector. Construct it with keyword
arguments; fields with a default value can be left out.

# Fields

- `id::Int64`: unique particle number. Ids must run from `1` to `N`, because each particle's
    random number generator is looked up by id.
- `type::Int64`: particle type. Forces and DOF evolvers only act on the types listed in
    their `ontypes` field.
- `m::Float64`: mass (**Default**: `1.0`)
- `zeta::Float64`: translational friction coefficient (**Default**: `1.0`)
- `zeta_R::Float64`: rotational friction coefficient (**Default**: `1.0`)
- `R::Float64`: radius (**Default**: `1.0`)
- `x::SVector{3,Float64}`: position, kept inside the box by periodic boundary conditions
- `xuw::SVector{3,Float64}`: unwrapped position, not affected by periodic boundary conditions.
    It is set equal to `x` at the start of a simulation. (**Default**: zeros)
- `v::SVector{3,Float64}`: velocity (**Default**: zeros)
- `f::SVector{3,Float64}`: total force. The forces add to it and the DOF evolvers reset it
    every timestep. (**Default**: zeros)
- `p::SVector{3,Float64}`: polarity vector (unit vector)
- `q::SVector{3,Float64}`: angular velocity; for particles moving in the xy-plane only the
    z-component is nonzero (**Default**: zeros)
- `T::SVector{3,Float64}`: total torque, accumulated and reset like `f` (**Default**: zeros)

# Example

```julia
p = Particles.Polar(id=1, type=1, R=1.0, x=[0.0, 0.0, 0.0], p=[1.0, 0.0, 0.0])
```
"""
@kwdef struct Polar <: Particle

    id::Int64
    type::Int64
    
    m::Float64   = 1.0
    zeta::Float64 = 1.0
    
    zeta_R::Float64 = 1.0

    R::Float64   = 1.0

    x::SVector{3, Float64}
    xuw::SVector{3, Float64} = @SVector [0.0, 0.0, 0.0]

    v::SVector{3, Float64}   = @SVector [0.0, 0.0, 0.0]
    f::SVector{3, Float64}   = @SVector [0.0, 0.0, 0.0]

    p::SVector{3, Float64}                                  #polarity vector
    q::SVector{3, Float64}   = @SVector [0.0, 0.0, 0.0]     #angular velocity (so for particles in xy this is in z)
    T::SVector{3, Float64}   = @SVector [0.0, 0.0, 0.0]     #torque  (so for particles in xy this is in z)
end

#Parametrize on the kernel length
@kwdef struct PolarKernel{Nk} <: Particle

    id::Int64
    type::Int64
    
    m::Float64   = 1.0
    zeta::Float64 = 1.0
    
    zeta_R::Float64 = 1.0

    R::Float64   = 1.0

    x::SVector{3, Float64}
    xuw::SVector{3, Float64} = @SVector [0.0, 0.0, 0.0]

    v::SVector{3, Float64}   = @SVector [0.0, 0.0, 0.0]
    f::SVector{3, Float64}   = @SVector [0.0, 0.0, 0.0]

    p::SVector{3, Float64}                                  #polarity vector
    q::SVector{3, Float64}   = @SVector [0.0, 0.0, 0.0]     #angular velocity (so for particles in xy this is in z)
    T::SVector{3, Float64}   = @SVector [0.0, 0.0, 0.0]     #torque  (so for particles in xy this is in z)

    kernel_x::SVector{Nk, Float64} = @SVector [0.0]
    kernel_y::SVector{Nk, Float64}= @SVector [0.0]

end
#Allow default kernel to work
PolarKernel(id, type, m, zeta, zeta_R, R, x, xuw, v, f, p, q, T, kernel_x, kernel_y) = PolarKernel{length(kernel_x)}(id, type, m, zeta, zeta_R, R, x, xuw, v, f, p, q, T, kernel_x, kernel_y)




"""
    PolarPolymer(; id, type, pol_id, id_in_pol, pol_N, x, p, kwargs...)

Polar particle that is one monomer of a polymer. It has all fields of [`Polar`](@ref),
plus three fields that describe its place in the polymer:

- `pol_id::Int64`: id of the polymer this particle belongs to
- `id_in_pol::Int64`: position of the particle along the polymer
- `pol_N::Int64`: number of particles in the polymer

`LPlot.polymers!` colours particles by polymer.
"""
@kwdef struct PolarPolymer <: Particle

    id::Int64
    type::Int64

    #id of the polymer this particle is part of
    pol_id::Int64

    #id/order within polymer
    id_in_pol::Int64

    #Number of particles in the polymer
    pol_N::Int64

    m::Float64   = 1.0
    zeta::Float64 = 1.0
    
    zeta_R::Float64 = 1.0

    R::Float64   = 1.0

    x::SVector{3, Float64}
    xuw::SVector{3, Float64} = @SVector [0.0, 0.0, 0.0]

    v::SVector{3, Float64}   = @SVector [0.0, 0.0, 0.0]
    f::SVector{3, Float64}   = @SVector [0.0, 0.0, 0.0]

    p::SVector{3, Float64}                                  #polarity vector
    q::SVector{3, Float64}   = @SVector [0.0, 0.0, 0.0]     #angular velocity (so for particles in xy this is in z)
    T::SVector{3, Float64}   = @SVector [0.0, 0.0, 0.0]     #torque  (so for particles in xy this is in z)
end



#Thanks to Nina Cankilic, Chinmay Pabshettiwar and Silke Henkes
#[Reference to MSc thesis to be inserted]
@kwdef struct GBEllipse<:Particle


    id::Int64
    type::Int64
    
    m::Float64   = 1.0
    zeta::Float64 = 1.0
    
    zeta_R::Float64 = 1.0

    R::Float64   = 1.0

    x::SVector{3, Float64}
    xuw::SVector{3, Float64} = @SVector [0.0, 0.0, 0.0]

    v::SVector{3, Float64}   = @SVector [0.0, 0.0, 0.0]
    f::SVector{3, Float64}   = @SVector [0.0, 0.0, 0.0]

    p::SVector{3, Float64}                                  #polarity vector
    q::SVector{3, Float64}   = @SVector [0.0, 0.0, 0.0]     #angular velocity (so for particles in xy this is in z)
    T::SVector{3, Float64}   = @SVector [0.0, 0.0, 0.0]     #torque  (so for particles in xy this is in z)

    stress::SMatrix{2,2,Float64,4} = @SMatrix [0.0 0.0 ; 0.0 0.0]

    Lambda::SMatrix{2,2,Float64,4}
    Lambda0::SMatrix{2,2,Float64,4}

    semi_major::Float64

    semi_minor::Float64

    semi_major0::Float64

    semi_minor0::Float64

    mu_alpha::Float64
    sigma_0::Float64 = 2.0

    aspect_ratio0::Float64

    K::Float64
    mu::Float64
    tau::Float64

end


#Parametrize on the kernel length and number of subparticles
@kwdef struct PolarShapeKernel{Nk, Ns} <: RigidBody

    id::Int64
    type::Int64
    
    m::Float64   = 1.0
    zeta::Float64 = 1.0
    
    zeta_R::Float64 = 1.0

    R::Float64   = 1.0

    x::SVector{3, Float64} #This now takes the interpretation as the COM location
    xuw::SVector{3, Float64} = @SVector [0.0, 0.0, 0.0]

    v::SVector{3, Float64}   = @SVector [0.0, 0.0, 0.0]
    f::SVector{3, Float64}   = @SVector [0.0, 0.0, 0.0]

    p::SVector{3, Float64}                                  #polarity vector
    q::SVector{3, Float64}   = @SVector [0.0, 0.0, 0.0]     #angular velocity (so for particles in xy this is in z)
    T::SVector{3, Float64}   = @SVector [0.0, 0.0, 0.0]     #torque  (so for particles in xy this is in z)

    kernel_x::SVector{Nk, Float64} = @SVector [0.0]
    kernel_y::SVector{Nk, Float64}= @SVector [0.0]

    # Initial positions of subparticles relative to COM, so will stay fixed because it is rigid
    xsrel0::SVector{Ns, Float64}
    ysrel0::SVector{Ns, Float64}

    #Updated relative positions of subparticles, so moves with p: absolute xs = xsrel + x
    xsrel::SVector{Ns, Float64} = xsrel0 .* p[1] .- ysrel0 .* p[2]
    ysrel::SVector{Ns, Float64} = xsrel0 .* p[2] .+ ysrel0 .* p[1]

    # Radii of subparticles
    Rs::SVector{Ns, Float64}

end

#Allow default kernel to work
PolarShapeKernel(id, type, m, zeta, zeta_R, R, x, xuw, v, f, p, q, T, kernel_x, kernel_y, xsrel0, ysrel0, xsrel, ysrel, Rs) = PolarShapeKernel{length(kernel_x),length(xsrel0)}(id, type, m, zeta, zeta_R, R, x, xuw, v, f, p, q, T, kernel_x, kernel_y, xsrel0, ysrel0, xsrel, ysrel, Rs)


#module end
end



