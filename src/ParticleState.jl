
# For single particle type, we can use StructArrays, which should eventually allow a GPU port. For multiple particle type mix, we need to fallback to union spltting. 
# We define the particle state function to handle this automatically for the user.

export ParticleState

"""
    ParticleState(particles)

Turn a vector of particles into the state a [`System`](@ref) needs.

- One particle type: returns a `StructArray` (one column per field, fastest).
- Several particle types: returns a `Vector{Union{T1, T2, ...}}`. Forces and DOF evolvers
  work unchanged, but the force loops are a few times slower.

Particle ids must run from `1` to `N` over the whole state.
"""
function ParticleState(particles::AbstractVector)
    isempty(particles) && error("JAMS: ParticleState needs at least one particle.")
    Ts = unique(typeof.(particles))
    length(Ts) == 1 && return StructArray(convert(Vector{only(Ts)}, particles))
    length(Ts) > 4 && @warn "JAMS: mixing $(length(Ts)) particle types. Keep in mind that Julia stops union-splitting for more than 4 types and the simulation becomes much slower."
    return convert(Vector{Union{Ts...}}, particles)
end

ParticleState(state::StructArray) = state

struct PRow{S}
    sa::S
    i::Int
end

# One particle type (StructArray): read and write single columns
@inline Base.getproperty(r::PRow{<:StructArray}, s::Symbol) =
    @inbounds getproperty(getfield(r, :sa), s)[getfield(r, :i)]
@inline function Base.setproperty!(r::PRow{<:StructArray}, s::Symbol, v)
    @inbounds getproperty(getfield(r, :sa), s)[getfield(r, :i)] = v
    return v
end

# Mixed types (Vector of immutable particles): write back a copy with one field changed
@inline Base.getproperty(r::PRow{<:AbstractVector}, s::Symbol) =
    @inbounds getproperty(getfield(r, :sa)[getfield(r, :i)], s)
@inline function Base.setproperty!(r::PRow{<:AbstractVector}, s::Symbol, v)
    sa = getfield(r, :sa)
    i = getfield(r, :i)
    @inbounds sa[i] = with_field(sa[i], Val(s), v)
    return v
end

#StructArray comes with its internal ability to set (change) immutable fields of structs. For the fallback option, we need to something similar: return a new struct with the field updated.
#It is generated, so it acts as if we wrote this manually for every particle type separately
@generated function with_field(p::T, ::Val{s}, v) where {T,s}
    s in fieldnames(T) || return :(error("JAMS: particle type ", $T, " has no field ", $(QuoteNode(s))))

    # Make new instance of the particle type, with the same fields as the old, except with the one new updated field
    args = [f === s ? :v : :(getfield(p, $(QuoteNode(f)))) for f in fieldnames(T)]

    # We are returning the generated function return, i.e. the return should be replaced by the actual particle type constructor.
    return :($T($(args...)))
end