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
struct Cells
    #For convenience store system sizes here as well
    sizes::NTuple{3,Float64}
    nbins::NTuple{3,Int} 
    lbins::NTuple{3,Float64}
    neighbors::Matrix{Int}
    # each column [:,i] corresponds to the cell ids of neighbours in the 27 directions
    # of cell i, a zero entry means no neighbour in that direction (in case of finite systems)
    # the entries are ordered so that every entry down the column after the first zero entry is also zero
    
    grouped_particle_ids::Vector{Int}

    group_edges::Vector{Int} #index in cell grouped_particle_ids at which a cell starts

    particle_id_to_cell_id::Vector{Int}

    counts::Vector{Int}         

    
end

@inline function find_cell_index(xi, nbins, lbins, sizes)
    cx = clamp(floor(Int, (xi[1] + sizes[1] / 2) / lbins[1]), 0, nbins[1] - 1) + 1
    cy = clamp(floor(Int, (xi[2] + sizes[2] / 2) / lbins[2]), 0, nbins[2] - 1) + 1
    cz = clamp(floor(Int, (xi[3] + sizes[3] / 2) / lbins[3]), 0, nbins[3] - 1) + 1
    return cx + nbins[1] * ((cy - 1) + nbins[2] * (cz - 1))

end

function construct_cell_neighbour_list(nbins, Periodic)

    cell_neighbour_list = zeros(Int, 27, prod(nbins))
    Nx = nbins[1]
    Ny = nbins[2]
    Nz = nbins[3]
    for cx in 1:Nx
        for cy in 1:Ny
            for cz in 1:Nz
                neighbour_number=1
                for x_offset in -1:1
                    for y_offset in -1:1
                        for z_offset in -1:1
                            cx_candidate = cx + x_offset
                            cy_candidate = cy + y_offset
                            cz_candidate = cz + z_offset

                            if Periodic
                                #mod 1 takes care of the index 1 convention in Julia 
                                # e.g. 10,10 stays at 10, while 11,10 becomes 1
                                cx_candidate = mod1(cx_candidate,Nx)
                                cy_candidate = mod1(cy_candidate,Ny)
                                cz_candidate = mod1(cz_candidate,Nz)

                            #Skip next part if not valid neigbour
                            elseif !( (1<=cx_candidate<=Nx) && (1<=cy_candidate<=Ny)  && (1<=cz_candidate<=Nz) )
                                continue
                            end
                            cell_candidate_index = cx_candidate + Nx* ( (cy_candidate-1) + Ny*(cz_candidate-1)) 
                            cell_index = cx + Nx*( (cy-1) + Ny*(cz-1))
                            #avoid duplicates in case of 1 single cell layer and pbc
                            @views if !(cell_candidate_index in cell_neighbour_list[:,cell_index])
                                cell_neighbour_list[neighbour_number,cell_index] = cell_candidate_index
                                neighbour_number+=1 #This way we fill in order and a 0 in the neighbour list then means no more neighbours
                            end

                        end
                    end
                end

            end
        end
    end
    return cell_neighbour_list
end


function construct_cell_list(system)
    #At least one bin
    rverlet = system.rcut_pair_global * (1 +system.skinfactor )
    nbins = ntuple(d -> max(floor(Int, system.sizes[d] / rverlet), 1), 3)
    lbins = ntuple(d -> system.sizes[d] / nbins[d], 3)

    N = length(system.initial_particle_state)
    Ncells = prod(nbins) #total number of bins
    
    #initialize cell list
    cells = Cells(system.sizes, nbins, lbins,construct_cell_neighbour_list(nbins, system.Periodic),
                  zeros(Int, N), zeros(Int, Ncells + 1), zeros(Int, N), zeros(Int, Ncells)
                  )

    return update_cells!(cells, system.initial_particle_state)
end

@inbounds function update_cells!(cells, current_particle_state)

    #Find new cell for each particle
    Threads.@threads for id in eachindex(current_particle_state)
        @inbounds cells.particle_id_to_cell_id[id] = find_cell_index(current_particle_state.x[id], cells.nbins, cells.lbins, cells.sizes)
    end

    #Reset counts
    fill!(cells.counts, 0)

    #Count how many particles in each cell
    for cell_id in cells.particle_id_to_cell_id
        cells.counts[cell_id] += 1
    end

    #Calculate cell-group edge indices based on counts
    cells.group_edges[1] = 1
    for i in eachindex(cells.counts)
         cells.group_edges[i+1] =  cells.group_edges[i] + cells.counts[i]
    end


    #Now we are going to place the particles id in the grouped-per-cell particle id list. We will reuse the counts array.
    #For each cell, the first particle id should go to the (left) group edge:
    for i in eachindex(cells.counts)
        cells.counts[i] = cells.group_edges[i]
    end

    for particle_id in eachindex(cells.particle_id_to_cell_id)

        #Particle id  needs to be placed in the group of cell id...
        cell_id = cells.particle_id_to_cell_id[particle_id]

        #Which for the first particle of the group is at
        cells.grouped_particle_ids[cells.counts[cell_id] ] = particle_id
        #and for later ones we need to shift to the right
        cells.counts[cell_id]+=1
 
    end


    return cells
end

struct Neighbours

    grouped_neighbour_ids::Vector{Int}

    group_edges::Vector{Int}

    counts::Vector{Int}  

    last_positions::Vector{SVector{3,Float64}}

    rverlet::Float64  

    rverlet2::Float64  

    skin::Float64
end

function construct_neighbour_list(cells, current_particle_state, system)


    Np = length(current_particle_state)
    
    skin = system.skinfactor * system.rcut_pair_global

    rverlet = system.rcut_pair_global + skin

    rverlet2 = rverlet^2

    neighbours = Neighbours(Int[], zeros(Int, Np + 1), zeros(Int, Np),zeros(SVector{3,Float64}, Np), rverlet, rverlet2, skin)
    
    return update_neighbour_list!(neighbours, cells, current_particle_state, system)

end

function update_neighbour_list!(neighbours, cells, current_particle_state, system)

    #First find length, then allocate.
    Threads.@threads for i in eachindex(current_particle_state)

        cell_id_of_p_i = cells.particle_id_to_cell_id[i]
        #Number of neighbours of particle i
        Nn=0
        @inbounds for k in 1:27
            candidate_cell_id = cells.neighbors[k, cell_id_of_p_i]
            candidate_cell_id == 0 && break #if zero: there will be no neigbouring cells anymore so we can stop completely 

            for ind in cells.group_edges[candidate_cell_id]:(cells.group_edges[candidate_cell_id+1]-1)

                j = cells.grouped_particle_ids[ind]

                j == i && continue
                dx = minimal_image_difference(current_particle_state.x[i], current_particle_state.x[j], system.sizes, system.Periodic)
                dxn2 = sum(abs2, dx)
                if  dxn2<= neighbours.rverlet2
                    Nn += 1
                end
            end
        end
        #Set number of neighbours of particle i
        neighbours.counts[i] = Nn
    end

    neighbours.group_edges[1] = 1
    for i in eachindex(neighbours.counts)
        neighbours.group_edges[i+1] = neighbours.group_edges[i] + neighbours.counts[i]
    end

    #Now we know the number of neighbours of each particle and where to place its neighbours's ids:
    #We might need to resize it.
    #Note neighbours.group_edges[end] = 1 +  sum(neighbours.counts)
    resize!(neighbours.grouped_neighbour_ids, neighbours.group_edges[end] - 1)

    #Now we have the information to fill grouped_neighbour_ids:
    Threads.@threads for i in eachindex(current_particle_state)

        cell_id_of_p_i = cells.particle_id_to_cell_id[i]
        #start fo group of neighbours of particle i
        start_ind=neighbours.group_edges[i]
        @inbounds for k in 1:27
            candidate_cell_id = cells.neighbors[k, cell_id_of_p_i]
            candidate_cell_id == 0 && break #if zero: there will be no neigbouring cells anymore so we can stop completely 

            for ind in cells.group_edges[candidate_cell_id]:(cells.group_edges[candidate_cell_id+1]-1)

                j = cells.grouped_particle_ids[ind]

                j == i && continue
                dx = minimal_image_difference(current_particle_state.x[i], current_particle_state.x[j], system.sizes, system.Periodic)
                dxn2 = sum(abs2, dx)
                if  dxn2<= neighbours.rverlet2
                    neighbours.grouped_neighbour_ids[start_ind] = j
                    start_ind += 1
                end
            end
        end
    end
    copyto!(neighbours.last_positions, current_particle_state.x)
    return neighbours

end

function check_for_rebuild(neighbours, current_particle_state, system)
    max_displacement2 = (neighbours.skin / 2)^2
    #Loop over all position components
    rebuild = false
    @inbounds for i in eachindex(current_particle_state.x)
        dx = minimal_image_difference(neighbours.last_positions[i], current_particle_state.x[i], system.sizes, system.Periodic)
        dxn2 = sum(abs2, dx)
        if dxn2>max_displacement2
            rebuild = true
        end
    end
    return rebuild
end

