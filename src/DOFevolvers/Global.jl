
"""
    polymer_p_set(ontypes)

Global DOF evolver that sets the polarity `p` of every monomer along the local tangent of
its polymer. Each timestep it adds the unit vectors along the bonds to the next monomer
(`x_(n+1) - x_n`) and from the previous monomer (`x_n - x_(n-1)`) to `p`, and then
normalizes `p`. Because this is repeated every timestep, `p` follows the polymer shape and
points from monomer `1` towards monomer `pol_N`. Use it instead of an orientation evolver
such as [`overdamped_pqT_xyc`](@ref); it does not use or reset the torque `T`.

# Fields

- `ontypes`: particle type (`Int`) or types (`Vector{Int}`) this evolver acts on. Both
    monomers of a bond must have a type in `ontypes`.

Requires particle fields `type`, `pol_id`, `id_in_pol`, `x` and `p`, e.g.
[`Particles.PolarPolymer`](@ref JAMS.Particles.PolarPolymer).
"""
struct polymer_p_set<: GlobalDOFevolver
    ontypes::Union{Int64,Vector{Int64}}
end
struct PRow{S}
    sa::S
    i::Int
end
@inline Base.getproperty(r::PRow, s::Symbol) =
    @inbounds getproperty(getfield(r, :sa), s)[getfield(r, :i)]
@inline function Base.setproperty!(r::PRow, s::Symbol, v)
    @inbounds getproperty(getfield(r, :sa), s)[getfield(r, :i)] = v
    return v
end


function evolve_globally!(current_particle_state, current_field_state, system, neighbours, dt, dofevolver::polymer_p_set)

    Threads.@threads for i in eachindex(current_particle_state)

        p_i =  PRow(current_particle_state, i)
        #Reset
        p_i.p = p_i.p*0

        @inbounds for n in neighbours.group_edges[i]:(neighbours.group_edges[i+1]-1)

            j = neighbours.grouped_neighbour_ids[n]
            
            p_j =  current_particle_state[j]

                    
            if p_i.type in dofevolver.ontypes && p_j.type in dofevolver.ontypes

                if p_i.pol_id==p_j.pol_id


                    if p_j.id_in_pol==p_i.id_in_pol+1
                        dx = minimal_image_difference(p_i.x, p_j.x, system.sizes, system.Periodic)

                        dxn = norm(dx)

                        p_i.p+=dx/dxn

                    elseif  p_j.id_in_pol==p_i.id_in_pol-1

                        dx = minimal_image_difference(p_i.x, p_j.x, system.sizes, system.Periodic)

                        dxn = norm(dx)

                        p_i.p+= -dx/dxn
                        end
                    end
                end
            end
        p_i.p= normalize(p_i.p)
    end
    return current_particle_state, current_field_state
end