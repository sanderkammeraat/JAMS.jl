module DOFevolvers
using StaticArrays
using LinearAlgebra
using Distributions
using Random
using SparseArrays


"""
    DOFevolver

Abstract supertype of all degree-of-freedom (DOF) evolvers. A DOF evolver turns the forces
and torques accumulated during a timestep into new positions, velocities and orientations,
and then resets the forces and torques to zero. Put DOF evolver instances in the
`dofevolvers` tuple of a [`System`](@ref JAMS.System).
"""
abstract type DOFevolver end

"""
    LocalDOFevolver <: DOFevolver

Supertype of DOF evolvers that update one particle at a time, using only that particle's
own fields. A new local DOF evolver needs a method

```julia
evolve_locally!(p_i, t, dt, dofevolver::MyEvolver)
```

that updates the fields of `p_i`, resets the forces and torques it used, and returns `p_i`.
"""
abstract type LocalDOFevolver <: DOFevolver end
include("DOFevolvers/Local.jl")

"""
    GlobalDOFevolver <: DOFevolver

Supertype of DOF evolvers that need the state of all particles at once. 
"""
abstract type GlobalDOFevolver <: DOFevolver end
#include("DOFevolvers/Global.jl")

"""
    FieldDOFevolver <: DOFevolver

Supertype of DOF evolvers for fields.
"""
abstract type FieldDOFevolver <: DOFevolver end
#include("DOFevolvers/Field.jl")


end