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

@inbounds function minimal_image_difference(xi, xj, system_sizes, system_Periodic)
    
    dx_x = minimal_image_difference_component(xj[1]-xi[1],system_sizes[1], system_Periodic)
    dx_y = minimal_image_difference_component(xj[2]-xi[2],system_sizes[2], system_Periodic)
    dx_z = minimal_image_difference_component(xj[3]-xi[3],system_sizes[3], system_Periodic)
    return SVector{3, Float64}(dx_x, dx_y, dx_z)
end

function minimal_image_difference_component(linear_difference,linear_size, system_Periodic)

        if system_Periodic
            if linear_difference>linear_size/2
                linear_difference-=linear_size
            end
            if linear_difference<=-linear_size/2
                linear_difference+=linear_size
            end
        end 

    return linear_difference
end


function evolve_globally!(current_particle_state, current_field_state, system, cells, dt, dofevolver::polymer_p_set)

    Threads.@threads for i in eachindex(current_particle_state)

        p_i =  PRow(current_particle_state, i)

        cell_id_of_p_i = cells.particle_id_to_cell_id[i]

        @inbounds for k in 1:27
            candidate_cell_id = cells.neighbors[k, cell_id_of_p_i]
            candidate_cell_id == 0 && break #if zero: there will be no neigbouring cells anymore so we can stop completely 

            for ind in cells.group_edges[candidate_cell_id]:(cells.group_edges[candidate_cell_id+1]-1)

                j = cells.grouped_particle_ids[ind]

                j == i && continue #If the particle in the candidate cell is particle i, skip the next part

                p_j =  PRow(current_particle_state, j)
                    
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
        end
    p_i.p= normalize(p_i.p)
    end
    return current_particle_state, current_field_state
end