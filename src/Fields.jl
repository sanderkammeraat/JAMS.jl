module Fields
using StaticArrays
using LinearAlgebra
using Distributions
using Random
using SparseArrays


abstract type Field end

struct ScalarField2d{Tbc}<:Field
    id::Int64
    type::Int64
    bin_centers::Tbc
    C::Matrix{Float64}
    Cv::Matrix{Float64}
    Cf::Matrix{Float64}
    lbin::Float64
end




#module end
end