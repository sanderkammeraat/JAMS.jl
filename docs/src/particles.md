# Particles

A particle type is an immutable struct that holds the degrees of freedom and parameters
of one particle. All particles of a simulation are stored together in a
[`ParticleState`](@ref), with one column per field, so `state.x` is the vector of all
positions.

Create particles with keyword arguments and collect them in a `ParticleState`:

```julia
initial_state = ParticleState([
    Particles.Polar(id=i, type=1, x=[0.0, 0.0, 0.0], p=[1.0, 0.0, 0.0])
    for i in 1:N
])
```

Particle ids must run from `1` to `N`. Vectors such as `x` always have three components,
also in 2D simulations.

To define your own particle type, see [A custom particle type](@ref).

## Available particles

| Particle | Description |
|---|---|
| [`Particles.Polar`](@ref) | Spherical particle with a position and a polarity vector. The standard choice for active Brownian particles and soft disks. |
| [`Particles.PolarPolymer`](@ref) | `Polar` particle that is one monomer of a polymer, with the polymer id, its position along the polymer and the polymer length. |

## Abstract types

| Type | Description |
|---|---|
| [`Particles.Particle`](@ref) | Supertype of all particle types. |
| [`Particles.RigidBody`](@ref) | Supertype of extended (rigid-body) particles. |

## Reference

```@docs
Particles.Polar
Particles.PolarPolymer
Particles.Particle
Particles.RigidBody
```
