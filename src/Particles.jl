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

Abstract supertype for extended (rigid-body) particles. Periodic boundary conditions are
only applied to their center position `x`.
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



    

#module end
end



