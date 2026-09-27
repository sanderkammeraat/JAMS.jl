module Forces
using StaticArrays
using LinearAlgebra
using Distributions
using Random
using SparseArrays
"""
    Force

Abstract supertype of all forces. Put force instances in the `forces` tuple of a
[`System`](@ref JAMS.System). The engine sorts them by subtype ([`ExternalForce`](@ref),
[`PairForce`](@ref) or [`FieldForce`](@ref)) and calls the matching method every timestep.
"""
abstract type Force end

"""
    ExternalForce <: Force

Supertype of forces and torques that act on one particle at a time, such as
self-propulsion or rotational noise. A new external force needs a method

```julia
contribute_external_force!(p_i, t, dt, rngs_particles, system, force::MyForce)
```

that adds to `p_i.f` and/or `p_i.T` and returns `p_i`. Here `p_i` is a mutable view of one
particle, `t` the current time, `dt` the timestep, `rngs_particles[p_i.id]` the random number
generator of this particle and `system` the [`System`](@ref JAMS.System).
"""
abstract type ExternalForce <:Force end
include("Forces/External.jl")

"""
    PairForce <: Force

Supertype of forces between pairs of particles. A new pair force needs a method

```julia
contribute_pair_force!(p_i, p_j, dx, dxn, t, dt, rngs_particles, system, force::MyForce)
```

that adds the force of particle `j` on particle `i` to `p_i.f` (and optionally a torque to
`p_i.T`) and returns `p_i`. Here `dx = x_j - x_i` is the separation vector (minimal image
if the system is periodic) and `dxn` its length. The method is only called for pairs with
`dxn <= system.rcut_pair_global`. It is called separately for `(i, j)` and `(j, i)`, so it
must only change `p_i`.
"""
abstract type PairForce <:Force end
include("Forces/Pair.jl")

"""
    FieldForce <: Force

Supertype of forces between particles and fields.
"""
abstract type FieldForce <:Force end
#include("Forces/Field.jl")

end