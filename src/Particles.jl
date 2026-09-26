module Particles

using StaticArrays
using LinearAlgebra
using Distributions
using Random
using SparseArrays

abstract type Particle end

abstract type RigidBody end

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



