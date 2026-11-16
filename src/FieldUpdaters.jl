module FieldUpdaters
using StaticArrays
using LinearAlgebra
using Distributions
using Random
using SparseArrays


# Fields have ghost cells. These ghost cells allow setting periodic or fixed boundary conditions. It also makes operations homogenous across inner field points.


abstract type FieldUpdater end

struct Diffusion<:FieldUpdater
    ontypes::Union{Int64,Vector{Int64}}
    D::Float64
end

function contribute_field_update!(field_i, t, dt, field_updater::Diffusion, rngs_fields)
 
    if field_i.type in field_updater.ontypes

        Dfactor = field_updater.D/field_i.lbin^2

        C = field_i.C
        Cf = field_i.Cf


        Threads.@threads for j in 2:size(field_i.C,2)-1

            @inbounds for i in 2:size(field_i.C,1)-1

                Cf[i,j] +=  Dfactor * (C[i+1,j] + C[i-1,j]  + C[i,j+1] + C[i,j-1]  - 4 * C[i,j] )
            end
        end
    end
    return field_i
end





# Module end
end
