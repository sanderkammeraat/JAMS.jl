module DOFevolvers
using StaticArrays
using LinearAlgebra
using Distributions
using Random
using SparseArrays


abstract type DOFevolver end

abstract type LocalDOFevolver <: DOFevolver end
include("DOFevolvers/Local.jl")

abstract type GlobalDOFevolver <: DOFevolver end
#include("DOFevolvers/Global.jl")

abstract type FieldDOFevolver <: DOFevolver end
#include("DOFevolvers/Field.jl")


end