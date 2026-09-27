
"""
    repulsive_soft_disk(ontypes, karray)
    repulsive_soft_disk(; ontypes, karray)

Harmonic repulsion between overlapping particles. Two particles at distance `r` interact
when `r < R_i + R_j`, with potential energy `U = k/2 * (R_i + R_j - r)^2`, where
`k = karray[type_i, type_j]`. There is no force between particles that do not overlap.

# Fields

- `ontypes`: particle type (`Int`) or types (`Vector{Int}`) this force acts on. Both
    particles must have a type in `ontypes`.
- `karray`: stiffness, a matrix indexed by the two particle types, e.g. `[1.0 2.0; 2.0 2.0]`
    for types 1 and 2. A plain number only works if all particles have type `1`.

Requires particle fields `type`, `R` and `f`. Set `rcut_pair_global` in the
[`System`](@ref JAMS.System) to at least twice the largest radius.
"""
@kwdef struct repulsive_soft_disk{T1} <: PairForce
    ontypes::Union{Int64,Vector{Int64}}
    karray::T1
end

function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, force::repulsive_soft_disk)


    if p_i.type in force.ontypes && p_j.type in force.ontypes
        d2R = p_i.R+p_j.R
        if dxn < d2R

            p_i.f += force.karray[p_i.type,p_j.type] .* (dxn-d2R) .* dx/dxn
        end
    end
    return p_i

end


"""
    morse(ontypes, Dearray, aarray)
    morse(; ontypes, Dearray, aarray)

Morse potential between particles, attractive beyond contact and repulsive below it:
`U = De * (1 - exp(-a * (r - r_e)))^2`, with equilibrium distance `r_e = R_i + R_j`,
`De = Dearray[type_i, type_j]` and `a = aarray[type_i, type_j]`.

The force acts on every pair closer than `rcut_pair_global`, so that cutoff sets the range
of the attraction.

# Fields

- `ontypes`: particle type (`Int`) or types (`Vector{Int}`) this force acts on. Both
    particles must have a type in `ontypes`.
- `Dearray`: well depth `De`, a matrix indexed by the two particle types
- `aarray`: inverse width `a`, a matrix indexed by the two particle types

A plain number instead of a matrix only works if all particles have type `1`.

Requires particle fields `type`, `R` and `f`.
"""
@kwdef struct morse{T1, T2}<:PairForce
    ontypes::Union{Int64,Vector{Int64}}
    Dearray::T1
    aarray::T2
    
end


function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, force::morse)

    if p_i.type in force.ontypes && p_j.type in force.ontypes
        re = p_i.R+p_j.R

        a = force.aarray[p_i.type,p_j.type]
        De = force.Dearray[p_i.type,p_j.type]

        p_i.f += -2 * De*a*( exp(-2a*(dxn-re)) - exp(-a*(dxn-re)) ) * dx/dxn

    end
    return p_i

end




