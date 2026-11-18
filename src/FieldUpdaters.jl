module FieldUpdaters
using StaticArrays
using LinearAlgebra
using Distributions
using Random
using SparseArrays


# Fields have ghost cells. These ghost cells allow setting periodic or fixed boundary conditions. It also makes operations homogenous across inner field points.


abstract type FieldUpdater end

@kwdef struct Diffusion<:FieldUpdater
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

@kwdef struct GhostSetPeriodic<:FieldUpdater
    ontypes::Union{Int64,Vector{Int64}}
end

function contribute_field_update!(field_i, t, dt, field_updater::GhostSetPeriodic, rngs_fields)

    if field_i.type in field_updater.ontypes

        set_ghost_values!(field_i.C)

        set_ghost_values!(field_i.Cv)

        set_ghost_values!(field_i.Cf)

    end
    return field_i
end

function set_ghost_values!(C)

    #corner points and edges are correctly set, see ghost cell update comments in engine
    C[1,:] .= @view C[end-1,:]
    C[end,:] .= @view C[2,:]

    C[:,1] .= @view C[:,end-1]
    C[:,end] .= @view C[:,2]


    return C
end

@kwdef struct Relax<:FieldUpdater
    ontypes::Union{Int64,Vector{Int64}}
    Cset::Float64
    k::Float64

end

function contribute_field_update!(field_i, t, dt, field_updater::Relax, rngs_fields)

    if field_i.type in field_updater.ontypes

        C = field_i.C
        Cf = field_i.Cf

            Threads.@threads for j in 2:size(field_i.C,2)-1

                @inbounds for i in 2:size(field_i.C,1)-1

                    Cf[i,j] += field_updater.k * (field_updater.Cset - C[i,j])

                end

            end
    end
    return field_i
end



# Module end
end
