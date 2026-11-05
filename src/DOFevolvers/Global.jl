
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
@kwdef struct polymer_p_set<: GlobalDOFevolver
    ontypes::Union{Int64,Vector{Int64}}
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

@kwdef struct overdamped_pairdis{T1, T2}<: GlobalDOFevolver
    ontypes::Union{Int64,Vector{Int64}}
    etaarray::T1
    rfactarray::T2
    linearcutoff::Bool
end


function evolve_globally!(current_particle_state, current_field_state, system, neighbours, dt, dofevolver::overdamped_pairdis)

    #We are solving \zeta_i v_i + \sum_j w_{ij} [v_i - v_j ] = F_i
    # Where w_{ij} is a specific weight between particles i and j, it could e..g be a function based on their distance, or just a constant.
    # Since this equation couples all different particles, but the dependence on the velocity vectors is linear, we can write it as a matrix equation:
    # [\zeta I + K ] v = F. 
    #We combine the matrix and divide by the friction coefficients to get v directly: Mv = F/zeta_array , with  M = [ I + K/zeta_I ], 

    #There might be a way to reduce this to a N x N matrix, but for now we treat this fully expanded in components: the matrix equation has dimensions 3N x 3N.

    # Let's assume that in worst case, a particle interacts with 6 neighbours, then we have the following number of nonzero elements: 3N (from the diagonal) + N * 6 * 3 = 3N + 18N = 21N out of 9N^2 total elements
    # Hence the non-sparsity of the matrix is 7 1/N. For 1000 particles this means a non-sparsity of only 0.7 %. Hence we should treat this as a sparse problem.

    # As outlined here: https://sparsearrays.juliasparse.org/dev/#man-sparse-performance, we should first collect nonzero indices and their values and construct the matrix afterwards.

    N = length(current_particle_state)
    Mdim = 3N

    row_inds = sizehint!(Int64[],  21*N)
    col_inds = sizehint!(Int64[],  21*N)
    vals = sizehint!(Float64[],  21*N) 

    counter=0
    for i in eachindex(current_particle_state)

        p_i = current_particle_state[i]

        #Fill self part, need to always do this for unaffected particles too
        for k in 1:3
            push!(row_inds, 3*(i - 1) + k )

            push!(col_inds, 3*(i - 1) + k )

            push!(vals, 1.)
        end

        if p_i.type in dofevolver.ontypes

            counter+=1
            zeta_i = p_i.zeta 
            
            @inbounds for n in neighbours.group_edges[i]:(neighbours.group_edges[i+1]-1)

                j = neighbours.grouped_neighbour_ids[n]

                p_j =  current_particle_state[j]

                if p_j.type in dofevolver.ontypes

                    dx = minimal_image_difference(p_i.x, p_j.x, system.sizes, system.Periodic)

                    dxn2 = sum(abs2,dx)

                    d2a = p_i.R+p_j.R

                    rfact = dofevolver.rfactarray[get_param_ind(dofevolver.ontypes,p_i.type),get_param_ind(dofevolver.ontypes,p_j.type)]
                    
                    if dxn2<= rfact* d2a

                        eta = dofevolver.etaarray[get_param_ind(dofevolver.ontypes,p_i.type),get_param_ind(dofevolver.ontypes,p_j.type)]

                        #dxn = sqrt(dxn2)
                        #Change below for different weights, e.g. based on distance
                        weight = dofevolver.linearcutoff ? 1. - sqrt(dxn2)/(rfact * d2a) : 1.

                        # val is multiplying v_i, so negate it for v_j
                        val = weight * eta / zeta_i
                        

                        for k in 1:3

                            # Fill neighbour parts
                            push!(row_inds, 3*(i - 1) + k )

                            push!(col_inds, 3*(j - 1) + k )

                            push!(vals, -val)

                            #Fill neighbour-self part
                            push!(row_inds, 3*(i - 1) + k )

                            push!(col_inds, 3*(i - 1) + k )

                            push!(vals, val)

                        end

                    end
                end
            end
        end
    end
    
    if counter>0
        # We have now prepared the values for M, let's construct it
        M = sparse(row_inds, col_inds, vals, Mdim, Mdim)

        sol = M \ vec(stack(current_particle_state.f))

        Threads.@threads for i in eachindex(current_particle_state)

            p_i =  PRow(current_particle_state, i)
            if p_i.type in dofevolver.ontypes
                p_i.v = @views sol[ 3(i-1)+1:3i]
                
                p_i.x += p_i.v * dt
                p_i.xuw += p_i.v * dt
                
                #reinitialize
                p_i.f *=0.
            end
        end
    end 
    return current_particle_state, current_field_state


end


@kwdef struct overdamped_polymer_pairdis{T1, T2}<: GlobalDOFevolver
    ontypes::Union{Int64,Vector{Int64}}
    etaarray::T1
    rfactarray::T2
    linearcutoff::Bool
end


function evolve_globally!(current_particle_state, current_field_state, system, neighbours, dt, dofevolver::overdamped_polymer_pairdis)

    #We are solving \zeta_i v_i + \sum_j w_{ij} [v_i - v_j ] = F_i
    # Where w_{ij} is a specific weight between particles i and j, it could e..g be a function based on their distance, or just a constant.
    # Since this equation couples all different particles, but the dependence on the velocity vectors is linear, we can write it as a matrix equation:
    # [\zeta I + K ] v = F. 
    #We combine the matrix and divide by the friction coefficients to get v directly: Mv = F/zeta_array , with  M = [ I + K/zeta_I ], 

    #There might be a way to reduce this to a N x N matrix, but for now we treat this fully expanded in components: the matrix equation has dimensions 3N x 3N.

    # Let's assume that in worst case, a particle interacts with 6 neighbours, then we have the following number of nonzero elements: 3N (from the diagonal) + N * 6 * 3 = 3N + 18N = 21N out of 9N^2 total elements
    # Hence the non-sparsity of the matrix is 7 1/N. For 1000 particles this means a non-sparsity of only 0.7 %. Hence we should treat this as a sparse problem.

    # As outlined here: https://sparsearrays.juliasparse.org/dev/#man-sparse-performance, we should first collect nonzero indices and their values and construct the matrix afterwards.

    N = length(current_particle_state)
    Mdim = 3N

    row_inds = sizehint!(Int64[],  21*N)
    col_inds = sizehint!(Int64[],  21*N)
    vals = sizehint!(Float64[],  21*N) 
    counter = 0

    for i in eachindex(current_particle_state)

        p_i = current_particle_state[i]

        #Fill self part, need to always do this for uncoupled particles too
        for k in 1:3
            push!(row_inds, 3*(i - 1) + k )

            push!(col_inds, 3*(i - 1) + k )

            push!(vals, 1.)
        end

        if p_i.type in dofevolver.ontypes

            counter+=1

            zeta_i = p_i.zeta 
            
            @inbounds for n in neighbours.group_edges[i]:(neighbours.group_edges[i+1]-1)

                j = neighbours.grouped_neighbour_ids[n]

                p_j =  current_particle_state[j]

                if p_j.type in dofevolver.ontypes && (p_i.pol_id!=p_j.pol_id)

                    dx = minimal_image_difference(p_i.x, p_j.x, system.sizes, system.Periodic)

                    dxn2 = sum(abs2,dx)

                    d2a = p_i.R+p_j.R

                    rfact = dofevolver.rfactarray[get_param_ind(dofevolver.ontypes,p_i.type),get_param_ind(dofevolver.ontypes,p_j.type)]
                    
                    if dxn2<= rfact* d2a

                        eta = dofevolver.etaarray[get_param_ind(dofevolver.ontypes,p_i.type),get_param_ind(dofevolver.ontypes,p_j.type)]

                        #dxn = sqrt(dxn2)
                        #Change below for different weights, e.g. based on distance
                        weight = dofevolver.linearcutoff ? 1. - sqrt(dxn2)/(rfact * d2a) : 1.

                        # val is multiplying v_i, so negate it for v_j
                        val = weight * eta / zeta_i
                        

                        for k in 1:3

                            # Fill neighbour parts
                            push!(row_inds, 3*(i - 1) + k )

                            push!(col_inds, 3*(j - 1) + k )

                            push!(vals, -val)

                            #Fill neighbour-self part
                            push!(row_inds, 3*(i - 1) + k )

                            push!(col_inds, 3*(i - 1) + k )

                            push!(vals, val)

                        end

                    end
                end
            end

        end
    end
    

    # We have now prepared the values for M, let's construct it
    if counter>0
        M = sparse(row_inds, col_inds, vals, Mdim, Mdim)

        sol = M \ vec(stack(current_particle_state.f))

        Threads.@threads for i in eachindex(current_particle_state)
            p_i =  PRow(current_particle_state, i)
            if p_i.type in dofevolver.ontypes
                

                p_i.v = @views sol[ 3(i-1)+1:3i]
                
                p_i.x += p_i.v * dt
                p_i.xuw += p_i.v * dt
                
                #reinitialize
                p_i.f *=0.
            end
        end
    end

    return current_particle_state, current_field_state


end