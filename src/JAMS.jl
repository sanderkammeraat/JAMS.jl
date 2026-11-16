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
include("ParticleState.jl")
include("Forces.jl")
include("Fields.jl")
include("FieldUpdaters.jl")
include("Save.jl")
include("LPlot.jl")
include("Initial.jl")
include("NeighbourFinding.jl")
include("DOFevolvers.jl")


#Make structs available to user
@reexport using .Particles
@reexport using .Forces
@reexport using .DOFevolvers
@reexport using .Fields
@reexport using .FieldUpdaters
@reexport using .LPlot
@reexport using .Initial
@reexport using .Save




include("Engine.jl")

end
