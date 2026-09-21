module Forces
using StaticArrays
using LinearAlgebra
using Distributions
using Random
using SparseArrays
abstract type Force end

abstract type ExternalForce <:Force end
include("Forces/External.jl")

abstract type PairForce <:Force end
include("Forces/Pair.jl")

abstract type FieldForce <:Force end
#include("Forces/Field.jl")

end