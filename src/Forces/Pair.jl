
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

"""
    polymer_harmonic_stretch(ontypes, internal_repulsion, karray, farray)

Harmonic bond between consecutive monomers of the same polymer (same `pol_id`,
`id_in_pol` differing by one): `U = k/2 * (r - l)^2`, with `k = karray[type_i, type_j]`.
The equilibrium bond length is `f * (R_i + R_j)`, with
`f = farray[type_i, type_j]`; choose `f < 1` for overlapping monomers.

If bonded monomers also repel each other through a soft-disk force of the same stiffness
(e.g. [`repulsive_soft_disk`](@ref)), set `internal_repulsion = true`. The rest length of
the spring is then shifted to `l = (2f - 1) * (R_i + R_j)`, so that spring and repulsion
together still settle at `f * (R_i + R_j)`. With [`polymer_repulsive_soft_disk`](@ref),
which skips bonded monomers, use `internal_repulsion = false` (then `l = f * (R_i + R_j)`).

# Fields

- `ontypes`: particle type (`Int`) or types (`Vector{Int}`) this force acts on. Both
    particles must have a type in `ontypes`.
- `internal_repulsion`: `true` if bonded monomers also feel a soft-disk repulsion (see above)
- `karray`: bond stiffness, a matrix indexed by the two particle types
- `farray`: bond length as a fraction of `R_i + R_j`, a matrix indexed by the two particle types

A plain number instead of a matrix only works if all particles have type `1`.

Requires particle fields `type`, `pol_id`, `id_in_pol`, `R` and `f`, e.g.
[`Particles.PolarPolymer`](@ref JAMS.Particles.PolarPolymer). Set `rcut_pair_global`
larger than the longest bond.
"""
struct polymer_harmonic_stretch{T1,T2}<:PairForce
    ontypes::Union{Int64,Vector{Int64}}
    internal_repulsion::Bool
    karray::T1
    farray::T2
end

function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, force::polymer_harmonic_stretch)

    if p_i.type in force.ontypes && p_j.type in force.ontypes

        #If in same polymer
        if p_i.pol_id==p_j.pol_id

            #If neigbouring points in the polymer
            if p_j.id_in_pol==p_i.id_in_pol+1 || p_j.id_in_pol==p_i.id_in_pol-1

                d2R = p_i.R+p_j.R

                f_factor = force.farray[p_i.type,p_j.type]
                
                if force.internal_repulsion
                    l_stretch = d2R*(2*f_factor - 1)
                else
                    l_stretch = d2R*f_factor
                end

                p_i.f+= force.karray[p_i.type,p_j.type] * (dxn-l_stretch) * dx/dxn
            end
        end
    end
    return p_i
end

"""
    polymer_harmonic_bend(ontypes, karray)

Bending stiffness of a polymer, from the discrete bending energy

`U = k/4 * Σ_n |R_{n+1} - R_{n}|^2`,

summed over the inner monomers `n = 0, …, pol_N - 2` of each polymer, with
`k = karray[type_i, type_j]`. For bonds of fixed length `b` this is
`U = k b^2 / 2 * Σ_n (1 - cos θ_n)`, with `θ_n` the angle between consecutive bonds, so a
straight polymer has the lowest energy. The force on a monomer comes from the monomers up
to two places away along the same polymer, including the special cases at the two ends.
Only use this for polymers with 5 or more monomers!

# Fields

- `ontypes`: particle type (`Int`) or types (`Vector{Int}`) this force acts on. Both
    particles must have a type in `ontypes`.
- `karray`: bending stiffness, a matrix indexed by the two particle types. A plain number
    only works if all particles have type `1`.

Requires particle fields `type`, `pol_id`, `id_in_pol`, `pol_N` and `f`, e.g.
[`Particles.PolarPolymer`](@ref JAMS.Particles.PolarPolymer), and polymers of at least 4
monomers. Pair forces are only evaluated within `rcut_pair_global`, so set it larger than
the distance between second neighbours along a polymer (about two bond lengths).
"""
struct polymer_harmonic_bend{T1}<:PairForce
    ontypes::Union{Int64,Vector{Int64}}
    karray::T1
end

@views function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, force::polymer_harmonic_bend)

    if p_i.type in force.ontypes && p_j.type in force.ontypes

        #If in same polymer
        if p_i.pol_id==p_j.pol_id

            #If i is one of the middle 
            if p_i.id_in_pol>2 && p_i.id_in_pol<p_i.pol_N-1


                if p_j.id_in_pol==p_i.id_in_pol+2 || p_j.id_in_pol==p_i.id_in_pol-2

                    p_i.f+= - force.karray[p_i.type,p_j.type]/2 * dx

                elseif p_j.id_in_pol==p_i.id_in_pol+1 || p_j.id_in_pol==p_i.id_in_pol-1

                    p_i.f+= 4 * force.karray[p_i.type,p_j.type]/2 * dx

                end

            #If i is the left most particle
            elseif p_i.id_in_pol==1

                if p_j.id_in_pol==3

                    p_i.f+= - force.karray[p_i.type,p_j.type]/2 * dx

                elseif p_j.id_in_pol==2

                    p_i.f+= 2 * force.karray[p_i.type,p_j.type]/2 * dx

                end
        
            #If i is the right most particle
            elseif p_i.id_in_pol==p_i.pol_N

                if p_j.id_in_pol==p_i.pol_N-2

                    p_i.f+= - force.karray[p_i.type,p_j.type]/2 * dx

                elseif p_j.id_in_pol==p_i.pol_N-1

                    p_i.f+= 2 * force.karray[p_i.type,p_j.type]/2 * dx

                end

            
            #If i is the second left most particle
            elseif p_i.id_in_pol==2 

                if p_j.id_in_pol==4

                    p_i.f+= - force.karray[p_i.type,p_j.type]/2 * dx

                elseif p_j.id_in_pol==3

                    p_i.f+= 4 *  force.karray[p_i.type,p_j.type]/2 * dx

                elseif p_j.id_in_pol==1

                    p_i.f+= 2 *  force.karray[p_i.type,p_j.type]/2 * dx

                end

            #If i is the second right most particle
            elseif p_i.id_in_pol==p_i.pol_N-1 #2   2

                if p_j.id_in_pol==p_i.pol_N-3 #0   4

                    p_i.f+= - force.karray[p_i.type,p_j.type]/2 * dx

                elseif p_j.id_in_pol==p_i.pol_N-2 #1  3

                    p_i.f+= 4 *  force.karray[p_i.type,p_j.type]/2 * dx

                elseif p_j.id_in_pol==p_i.pol_N #3   1

                    p_i.f+= 2 *  force.karray[p_i.type,p_j.type]/2 * dx

                end
            end
        
        end
    end
    return p_i
end

"""
    polymer_pair_polar_nematic(ontypes, bundles, rfact, v0)

Active pair force between monomers of *different* polymers that are closer than
`r_c = rfact * (R_i + R_j)`. The strength decreases linearly with distance,
`β = 1 - r / r_c`, and the force is directed along the polarities `p` (for polymers, the
local tangent set by [`DOFevolvers.polymer_p_set`](@ref JAMS.DOFevolvers.polymer_p_set)):

`f_i += β v0 (p_i - p_j)`

If `bundles` is `true`, the antiparallel rule `f_i += β v0 (p_i - p_j)` is used for all
pairs, so parallel monomers feel no force and form nonmotile bundles.

However, if `bundles` is `false`, for parallel polymers, we use f_i += β v0 (p_i + p_j) * sign_ij, with 
sign_ij = sign(p_i.pol_id - p_j.pol_id). 

# Fields

- `ontypes`: particle type (`Int`) or types (`Vector{Int}`) this force acts on. Both
    particles must have a type in `ontypes`.
- `bundles`: if `true`, use the antiparallel rule for all pairs (see above)
- `rfact`: interaction range, as a multiple of `R_i + R_j`
- `v0`: force strength 

Requires particle fields `type`, `pol_id`, `R`, `p` and `f`, e.g.
[`Particles.PolarPolymer`](@ref JAMS.Particles.PolarPolymer). Set `rcut_pair_global` to at
least `rfact` times the largest `R_i + R_j`.
"""
struct polymer_pair_polar_nematic<:PairForce
    ontypes::Union{Int64,Vector{Int64}}
    bundles::Bool
    rfact::Float64
    v0::Float64
end
function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, force::polymer_pair_polar_nematic)
    
    if p_i.type in force.ontypes && p_j.type in force.ontypes

        if p_j.pol_id!=p_i.pol_id

            d2a = p_i.R+p_j.R
            r = force.rfact*d2a::Float64

            if dxn < r

                β = 1 - dxn/r

                if force.bundles

                    p_i.f += β*force.v0[p_i.type] * (p_i.p .- p_j.p)

                else

                    if dot(p_i.p, p_j.p) > 0

                        sign_ij = sign(p_i.pol_id - p_j.pol_id)

                        p_i.f += β*sign_ij*force.v0[p_i.type] * (p_i.p .+ p_j.p)

                    else

                        p_i.f += β*force.v0[p_i.type] * (p_i.p .- p_j.p)

                    end

                end
            end

        end

    end
    return p_i

end

"""
    polymer_repulsive_soft_disk(ontypes, karray)

Same harmonic repulsion as [`repulsive_soft_disk`](@ref), `U = k/2 * (R_i + R_j - r_ij)^2`
for overlapping particles, except between bonded monomers (consecutive monomers of the
same polymer), which may overlap freely. Use it together with
[`polymer_harmonic_stretch`](@ref) with `internal_repulsion = false`.

# Fields

- `ontypes`: particle type (`Int`) or types (`Vector{Int}`) this force acts on. Both
    particles must have a type in `ontypes`.
- `karray`: stiffness, a matrix indexed by the two particle types. A plain number only
    works if all particles have type `1`.

Requires particle fields `type`, `pol_id`, `id_in_pol`, `R` and `f`, e.g.
[`Particles.PolarPolymer`](@ref JAMS.Particles.PolarPolymer).
"""
struct polymer_repulsive_soft_disk{T1}<:PairForce
    ontypes::Union{Int64,Vector{Int64}}
    karray::T1
end
function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, force::polymer_repulsive_soft_disk)

    if p_i.type in force.ontypes && p_j.type in force.ontypes
        if !(p_i.pol_id==p_j.pol_id) || ((p_i.pol_id==p_j.pol_id) && abs(p_i.id_in_pol-p_j.id_in_pol) > 1)
        d2R = p_i.R+p_j.R
            if dxn < d2R

                p_i.f += force.karray[p_i.type,p_j.type] * (dxn-d2R) * dx/dxn

            end
        end
    end
    return p_i

end

"""
    polymer_pairAN(ontypes, torque, traceless, intrapol, rfact, k_par, k_per, parray)

[`pairAN`](@ref) between monomers of *different* polymers only. See [`pairAN`](@ref) for
the force and the meaning of `torque`, `traceless`, `rfact`, `k_par`, `k_per` and `parray`.

The field `intrapol` is stored but not used yet: monomers of the same polymer never
interact through this force.

Requires particle fields `type`, `pol_id`, `R`, `p`, `f` (and `T` if `torque` is `true`),
e.g. [`Particles.PolarPolymer`](@ref JAMS.Particles.PolarPolymer).
"""
struct polymer_pairAN{T1,T2, T3}<:PairForce
    ontypes::Union{Int64,Vector{Int64}}
    torque::Bool
    traceless::Bool
    intrapol::Bool
    rfact::Float64
    k_par::T1
    k_per::T2
    parray::T3
end
function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, force::polymer_pairAN)
    
    if p_i.type in force.ontypes && p_j.type in force.ontypes

        if  p_j.pol_id!=p_i.pol_id

            d2a = p_i.R+p_j.R
            

            r = force.rfact*d2a::Float64

            if dxn < r

                z_hat = @SVector [0,0,1]

                f = @MVector zeros(length(dx))
                β = 1 - dxn/r
                sigma_i_dot_dx = @MVector zeros(length(dx))

                sigma_j_dot_dx = @MVector zeros(length(dx))

                sigma_i_dot_dx_perp = @MVector zeros(length(dx))

                sigma_j_dot_dx_perp = @MVector zeros(length(dx))
                    
                if force.traceless
                    sigma_i_dot_dx = force.parray[p_i.type] * ( dot(p_i.p,dx) .* p_i.p - 0.5 * dx)

                    sigma_j_dot_dx = force.parray[p_j.type] * ( dot(p_j.p, dx) .* p_j.p - 0.5 * dx)



                    sigma_i_dot_dx_perp = force.parray[p_i.type] *  (dot(p_i.p,cross(dx,z_hat)) .* p_i.p - 0.5 * cross(dx,z_hat)) 

                    sigma_j_dot_dx_perp = force.parray[p_j.type] *  (dot(p_j.p,cross(dx,z_hat)).* p_j.p- 0.5 * cross(dx,z_hat)) 


                else
                    #par
                    sigma_i_dot_dx = force.parray[p_i.type] * dot(p_i.p,dx) .* p_i.p 

                    sigma_j_dot_dx = force.parray[p_j.type] * dot(p_j.p, dx) .* p_j.p

                    # perp
                    sigma_i_dot_dx_perp = force.parray[p_i.type] *  dot(p_i.p,cross(dx,z_hat)) .* p_i.p

                    sigma_j_dot_dx_perp = force.parray[p_j.type] *  dot(p_j.p,cross(dx,z_hat)).* p_j.p

                end


                f = β .*  ( force.k_par .* (sigma_i_dot_dx .+ sigma_j_dot_dx )  .+ force.k_per .* (sigma_i_dot_dx_perp .+ sigma_j_dot_dx_perp ))

                p_i.f += f
                #add torque
                if force.torque
                    p_i.T+= cross(dx/2, f)
                end
            end
        end
    end
    return p_i

end

"""
    pairAN(ontypes, torque, traceless, rfact, k_par, k_per, parray)

Pairwise force from active nematic stresses, between particles closer than
`r_c = rfact * (R_i + R_j)`. Every particle carries an active stress
`σ = s_i Q_i`, with
strength `s_i = parray[type]`. The force on `i` is

`f_i += β [k_par (σ_i + σ_j) ⋅ dx + k_per (σ_i + σ_j) ⋅ (dx × ẑ)]`,

with `dx = x_j - x_i` and `β = 1 - r / r_c`.

If `torque` is `true`, the force also gives a torque `T_i += (dx / 2) × f_i`, as if it
acted at the midpoint between the particles.

The traceless stress and the perpendicular direction `dx × ẑ` assume particles in
the xy-plane.

# Fields

- `ontypes`: particle type (`Int`) or types (`Vector{Int}`) this force acts on. Both
    particles must have a type in `ontypes`.
- `torque`: if `true`, also add the torque `(dx / 2) × f_i`
- `traceless`: if `true`, use the traceless stress `s (p p - I/2)`
- `rfact`: interaction range, as a multiple of `R_i + R_j`
- `k_par`: strength of the force along `dx` (a number)
- `k_per`: strength of the force perpendicular to `dx` (a number)
- `parray`: active stress strength `p`, a vector indexed by particle type. 

Requires particle fields `type`, `R`, `p`, `f` (and `T` if `torque` is `true`). Set
`rcut_pair_global` to at least `rfact` times the largest `R_i + R_j`.
"""
struct pairAN{T1,T2, T3}<: PairForce
    ontypes::Union{Int64,Vector{Int64}}
    torque::Bool
    traceless::Bool
    rfact::Float64
    k_par::T1
    k_per::T2
    parray::T3
end



function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, force::pairAN)
    
    if p_i.type in force.ontypes && p_j.type in force.ontypes
        d2a = p_i.R+p_j.R

        r = force.rfact*d2a::Float64

        if dxn < r

            z_hat = @SVector [0,0,1]

            f = @MVector zeros(length(dx))

            sigma_i_dot_dx = @MVector zeros(length(dx))

            sigma_j_dot_dx = @MVector zeros(length(dx))

            sigma_i_dot_dx_perp = @MVector zeros(length(dx))

            sigma_j_dot_dx_perp = @MVector zeros(length(dx))

            β = 1 - dxn/r
            #rij_cap = dx/dxn


            if force.traceless
                sigma_i_dot_dx.= force.parray[p_i.type] .* ( dot(p_i.p,dx) .* p_i.p - 0.5 .* dx)

                sigma_j_dot_dx.= force.parray[p_j.type] .* ( dot(p_j.p, dx) .* p_j.p - 0.5 .* dx)



                sigma_i_dot_dx_perp.= force.parray[p_i.type] .*  (dot(p_i.p, cross(dx,z_hat)) .* p_i.p - 0.5 .* cross(dx, z_hat)) 

                sigma_j_dot_dx_perp.= force.parray[p_j.type] .*  (dot(p_j.p, cross(dx,z_hat)).* p_j.p- 0.5 .* cross(dx, z_hat)) 


            else
                #par
                sigma_i_dot_dx.= force.parray[p_i.type] .* dot(p_i.p,dx) .* p_i.p 

                sigma_j_dot_dx.= force.parray[p_j.type] .* dot(p_j.p, dx) .* p_j.p

                # perp
                sigma_i_dot_dx_perp.= force.parray[p_i.type] .*  dot(p_i.p,cross(dx,z_hat)) .* p_i.p

                sigma_j_dot_dx_perp.= force.parray[p_j.type] .*  dot(p_j.p,cross(dx,z_hat)).* p_j.p

            end
            f= β .*  ( force.k_par .* (sigma_i_dot_dx .+ sigma_j_dot_dx )  .+ force.k_per .* (sigma_i_dot_dx_perp .+ sigma_j_dot_dx_perp ))

            p_i.f+= f
            #add torque
            if force.torque
                p_i.T += cross(dx/2, f)
            end
        end
    end
    return p_i

end


"""
    pair_nematic_alignment(ontypes, rcut, J)

Nematic alignment torque between particles closer than `rcut`:
`T_i += J (p_i × p_j) (p_i ⋅ p_j)`. In the xy-plane this is `J/2 sin(2(θ_j - θ_i))`
about z, which aligns `p_i` with either `p_j` or `-p_j`, whichever is closer.

# Fields

- `ontypes`: particle type (`Int`) or types (`Vector{Int}`) this force acts on. Both
    particles must have a type in `ontypes`.
- `rcut`: interaction range. Pairs further apart than `rcut_pair_global` are never
    evaluated, so a larger `rcut` has no effect.
- `J`: alignment strength

Requires particle fields `type`, `p` and `T`.
"""
struct pair_nematic_alignment<: PairForce
    ontypes::Union{Int64,Vector{Int64}}
    rcut::Float64
    J::Float64
end
function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt, rngs_particles, system, force::pair_nematic_alignment)
    
    if p_i.type in force.ontypes && p_j.type in force.ontypes
        #d2a = p_i.R+p_j.R

        if dxn < force.rcut

            #add torque
            p_i.T+=force.J*cross(p_i.p, p_j.p) .* (dot(p_i.p,p_j.p))
        end
    end
    return p_i

end

"""
    pair_polar_alignment(ontypes, rcut, J)

Polar alignment torque between particles closer than `rcut`:
`T_i += J (p_i × p_j)`. In the xy-plane this is `J sin(θ_j - θ_i)` about z, which rotates
`p_i` towards `p_j`.

# Fields

- `ontypes`: particle type (`Int`) or types (`Vector{Int}`) this force acts on. Both
    particles must have a type in `ontypes`.
- `rcut`: interaction range. Pairs further apart than `rcut_pair_global` are never
    evaluated, so a larger `rcut` has no effect.
- `J`: alignment strength

Requires particle fields `type`, `p` and `T`.
"""
struct pair_polar_alignment<:PairForce
    ontypes::Union{Int64,Vector{Int64}}
    rcut::Float64
    J::Float64
end

function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt, rngs_particles, system, force::pair_polar_alignment)
    
    if p_i.type in force.ontypes && p_j.type in force.ontypes
        #d2a = p_i.R+p_j.R

        if dxn < force.rcut

            #add torque
            p_i.T+=force.J*cross(p_i.p, p_j.p)
        end
    end
    return p_i

end




