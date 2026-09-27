module JAMS

#Dependencies
using StaticArrays
using HDF5
using ProgressMeter
using SparseArrays
using Observables
using StructArrays
using Reexport
using CairoMakie
using JLD2
using Base.Cartesian

#Convenience export of packages for new users. Discutable practice...
@reexport using Random
@reexport using Distributions
@reexport using LinearAlgebra


#Stub function in case GLMakie is not used (e.g. on a headless cluster)
function setup_system_plotting end
function GLMakie_window_closeall end

include("Particles.jl")
include("Forces.jl")
include("DOFevolvers.jl")
include("Fields.jl")
include("FieldUpdaters.jl")
include("Save.jl")
include("LPlot.jl")
include("Initial.jl")


#Make structs available to user
@reexport using .Particles
@reexport using .Forces
@reexport using .DOFevolvers
@reexport using .Fields
@reexport using .FieldUpdaters
@reexport using .LPlot
@reexport using .Initial
@reexport using .Save

#Convenience alias for StructArray for new users
export ParticleState
"""
    ParticleState(particles)

Alias for `StructArray` from StructArrays.jl. Use it to turn a vector of particles into
the state that a [`System`](@ref) needs:

```julia
initial_state = ParticleState([Particles.Polar(id=i, type=1, x=..., p=...) for i in 1:N])
```

Each particle field becomes a column, so `initial_state.x` is the vector of all positions,
and `initial_state[i]` returns particle `i`.
"""
const ParticleState = StructArray


include("Engine.jl")

end
