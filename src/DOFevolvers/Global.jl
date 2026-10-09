
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