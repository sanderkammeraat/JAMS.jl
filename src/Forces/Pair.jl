
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




