@inline function periodic!(p_i, systemsizes)

    p_i.x = mod.(p_i.x .+ systemsizes ./ 2, systemsizes) .- systemsizes ./ 2
    return p_i
end

@inline function periodic!(p_i::Particles.RigidBody, systemsizes)

p_i.x = mod.(p_i.x .+ systemsizes ./ 2, systemsizes) .- systemsizes ./ 2
    #Deliberately not updating the extend points
    # for j=1:size(p_i.xe)[1]
    #     for (i, xi) in pairs(p_i.xe[j,:])

    #         if xi<-systemsizes[i]/2
    #             p_i.xe[j,i] = xi + systemsizes[i]
    #         elseif xi>=systemsizes[i]/2
    #             p_i.xe[j,i] = xi - systemsizes[i]
    #         end
    #     end
    
    # end    
    return p_i
end

#Initialize unwrapped coordinates to save the user the hassle to set equal to the initial wrapped coordinates
function init_unwrap!(p_i, t)

    if t==0
        p_i.xuw= copy(p_i.x)
    end

    return p_i
end
#Initialize forces and torques to zero, for easy chaining of sims
function init_f_T!(p_i, t)

    if t==0
        p_i.f*=0.
        p_i.T*=0.
    end
    return p_i
end


#Optimize for field indices
function minimal_image_closest_field_center(x, bin_centers, lbin)

    x_ind = length(bin_centers[1])>1 ? clamp(round(Int, (x[1] - bin_centers[1][2]) / lbin) + 2, 2, length(bin_centers[1])-1) : 1
    y_ind = length(bin_centers[2])>1 ? clamp(round(Int, (x[2] - bin_centers[2][2]) / lbin) + 2, 2, length(bin_centers[2])-1) : 1

    z_ind = length(bin_centers[3])>1 ? clamp(round(Int, (x[3] - bin_centers[3][2]) / lbin) + 2, 2, length(bin_centers[3])-1) : 1

    return  SVector{3, Int64}(x_ind, y_ind, z_ind)

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


export System
@kwdef struct System{Tips, Tifs, Tfor , Tfu , Tdof}

    #Vector that determines the linear size of the system
    sizes::NTuple{3, Float64}

    #Array containing particles in a specific state
    initial_particle_state::Tips

    #Array containing fields in a specific state
    initial_field_state::Tifs = ()

    #Array of force structs:
    forces::Tfor

    field_updaters::Tfu = ()
    
    #Array of structs to evolve dof (and reinitialize forces)
    dofevolvers::Tdof

    #Spatially periodic boundary conditions?
    Periodic::Bool

    #Global cutoff for pairwise interactions
    rcut_pair_global::Float64
end
#Output formatter to conveniently chain simulations in one .jl file without the need of intermediate saving to disk
struct SIM{T1, T2, T3, T4}
    final_particle_state::T1
    final_field_state::T2
    dt::T3
    t_stop::T4
    system::System #Note that system will contain the initial states
    finished::Bool
end


function save_raw_obj_data!(file_group, obj)

    obj_name = string(nameof(typeof(obj)))

    field_names = fieldnames(typeof(obj))

    obj_group = create_group(file_group, obj_name)

    for field_name in field_names

        name = string(field_name)
        val = getfield(obj, field_name)

        if val isa Bool
            obj_group[name] = string(val)
        elseif val isa StaticArray
            obj_group[name] = collect(val)
        else
            obj_group[name] = val
        end

    end

    return file_group
end


function save_raw_metadata!(file, system,external_forces,pair_forces, field_forces,local_dofevolvers, global_dofevolvers, field_dofevolvers,
     integration_tax,dt,t_stop,Tsave,save_tax, master_seed)

    create_group(file, "system")

    create_group(file["system"],"forces")

    create_group(file["system"]["forces"],"external")

    create_group(file["system"]["forces"],"pair")

    create_group(file["system"]["forces"],"field")


    for force in external_forces

        group = file["system"]["forces"]["external"]

        save_raw_obj_data!(group, force)

    end
    for force in pair_forces

        group = file["system"]["forces"]["pair"]

        save_raw_obj_data!(group, force)

    end
    for force in field_forces

        group = file["system"]["forces"]["field"]

        save_raw_obj_data!(group, force)

    end

    create_group(file["system"],"field_updaters")


    for fieldupdater in system.field_updaters

        group = file["system"]["field_updaters"]

        save_raw_obj_data!(group, fieldupdater)

    end

    create_group(file["system"],"dofevolvers")

    create_group(file["system"]["dofevolvers"],"local")

    create_group(file["system"]["dofevolvers"],"global")

    create_group(file["system"]["dofevolvers"],"field")


    for dofevolver in local_dofevolvers

        group = file["system"]["dofevolvers"]["local"]

        save_raw_obj_data!(group, dofevolver)
    end

    for dofevolver in global_dofevolvers

        group = file["system"]["dofevolvers"]["global"]

       save_raw_obj_data!(group, dofevolver)
    end

    for dofevolver in field_dofevolvers

        group = file["system"]["dofevolvers"]["field"]

        save_raw_obj_data!(group, dofevolver)
    end

    create_group(file, "integration_info")


    file["integration_info"]["integration_tax"] = integration_tax

    file["integration_info"]["Tsave"] = Tsave

    file["integration_info"]["save_tax"] = save_tax

    file["integration_info"]["dt"] = dt

    file["integration_info"]["t_stop"] = t_stop

    file["integration_info"]["master_seed"] = string(master_seed)


    file["system"]["sizes"] = collect(system.sizes)
    file["system"]["rcut_pair_global"] = system.rcut_pair_global
    file["system"]["periodic"] = string(system.Periodic)

    return file

end


export Euler_integrator 
"""
    Euler_integrator(system, dt, t_stop; kwargs...) -> SIM
Evolves the initial state defined in system according to the forces and DOFevolvers in it, using the Euler(-Maruyama) algorithm.
Returns SIM struct with final states for chaining simulations.

# Arguments

- `system`: Instance of System struct
- `dt`: Integration time step size
- `t_stop`: Integration end time

# Keywords

- `seed`: master seed for reproducibility
    (**Default**: `nothing`)
- `Tsave`: Save every Tsave dt-steps
    (**Default**: `nothing`)
- `save_functions`: Tuple of functions to save data
    (**Default**: `nothing`)
- `save_folder_path`: Folder path to save data (need not exist)
    (**Default**: `nothing`)
- `save_tag`: Add custom name to save file name, so that it becomes "save_tag.raw_data.h5"
    (**Default**: `nothing`)
- `Tplot`: Plot every Tplot dt-steps
    (**Default**: `nothing`)
- `fps`: Default frames per second for rendering of live plotting
    (**Default**: `30`)
- `plot_functions`: Tuple of plot function to live plot
    (**Default**: `nothing`)
- `plotdim`: Default plotting dimension, change to 3 for 3d plot
    (**Default**: `2`)
- `record_folder_path`: Folder path to save live recording in, need not exist
    (**Default**: `nothing`)
- `crf`: Default compression rate for saving the live recording
    (**Default**: `23`)
- `res`: Pixel resoluation of plotting, e.g. (1000,1000)
    (**Default**: `nothing`)
- `format`: Default format for saving live recording
    (**Default**: `"mp4"`)
- `sbs`: Set to true for side-by-side plotting, to be used i.c.m. with plotdim=3
    (**Default**: `false`)
"""
function Euler_integrator(system, dt, t_stop; seed=nothing, Tsave=nothing, save_functions=nothing, save_folder_path=nothing, save_tag=nothing, Tplot=nothing, fps=30, plot_functions=nothing,plotdim=2,record_folder_path=nothing,crf=23,res=nothing,format="mp4",sbs=false)


    integration_tax = collect(0:dt:t_stop)



    #This is really amazing: Julia rng is dependent on the task spawn structure, NOT on the 
    # parallel executation schedule, see https://julialang.org/blog/2021/11/julia-1.7-highlights/#new_rng_reproducible_rng_in_tasks
    # for more info, so this makes the program, even with mulitithreading and dynamic thread scheduling reproducible if we create a rng for every particle and field!
    # See for info https://github.com/JuliaLang/julia/issues/49064 !!
    if !isnothing(seed)
        master_seed = seed
        Random.seed!(seed)

    else
    #Otherwise generate one for reproducibility
    #Reset seed in case user used a specific seed for initial conditions
        Random.seed!()
        seed = rand(1:100000000000000000000000000000)
        master_seed = seed
        Random.seed!(seed)
    end
    



    if system.Periodic==false
        @warn ("JAMs: System is set to non-periodic: you should make sure particles always stay in system sizes for correctly working cell lists. The program will catch this by throwing an error if a particle is detected outside the box.")
    end

    cells = construct_cell_list(system)

    external_forces = Tuple(force for force in system.forces if isa(force, Forces.ExternalForce))
    pair_forces = Tuple(force for force in system.forces if isa(force, Forces.PairForce))

    field_forces = Tuple(force for force in system.forces if isa(force, Forces.FieldForce))


    Next = length(external_forces)
    Npair = length(pair_forces)

    Nfield = length(field_forces)

    Nfieldu = length(system.field_updaters)

    local_dofevolvers = Tuple(evolver for evolver in system.dofevolvers if isa(evolver, DOFevolvers.LocalDOFevolver))

    global_dofevolvers = Tuple(evolver for evolver in system.dofevolvers if isa(evolver, DOFevolvers.GlobalDOFevolver))

    field_dofevolvers = Tuple(evolver for evolver in system.dofevolvers if isa(evolver, DOFevolvers.FieldDOFevolver))

    

    current_particle_state = deepcopy(system.initial_particle_state)
    current_field_state = deepcopy(system.initial_field_state)



    rngs_fields = [Xoshiro(master_seed+i) for i in eachindex(current_field_state)]

    #Assuming number of fields stay constant
    rngs_particles = [Xoshiro(length(current_field_state)+master_seed+i) for i in eachindex(current_particle_state)]


    if !isnothing(Tsave)
        save_nax = [n for n in eachindex(integration_tax) if (n-1)%Tsave==0]
        save_tax = [ integration_tax[n] for n in eachindex(integration_tax) if (n-1)%Tsave==0 ]

        #Prepare save folder
        #If folder already exists, simply returns folder path. If folder is not existing, it will create the (sub)folders and return the path
        mkpath(save_folder_path)

        if isnothing(save_tag)
            JAMs_file_name = "JAMs_container.jld2"

            JAMs_final_state_file_name = "JAMs_final_state.jld2"
            
            raw_data_file_name = "raw_data.h5"
        else
            JAMs_file_name = save_tag * "_"* "JAMs_container.jld2"

            JAMs_final_state_file_name = save_tag * "_"* "JAMs_final_state.jld2"
            
            raw_data_file_name = save_tag * "_"*"raw_data.h5"
        end

        if isfile(joinpath(save_folder_path, JAMs_file_name)) || isfile(joinpath(save_folder_path, raw_data_file_name)) ||  isfile(joinpath(save_folder_path, JAMs_final_state_file_name))
            error("JAMs: Specified save folder already contains JAMs file(s): " * JAMs_file_name * " and/or " * raw_data_file_name*". JAMs aborted to prevent overwriting.")
        end

        #Store system and integration info in JAMS container
        jldopen(joinpath(save_folder_path, JAMs_file_name),"a+") do JAMs_file

            JAMs_file["system"] = system

            JAMs_file["integration_info/integration_tax"] = integration_tax

            JAMs_file["integration_info/dt"] = dt

            JAMs_file["integration_info/t_stop"] = t_stop

            JAMs_file["integration_info/Tsave"] = Tsave

            JAMs_file["integration_info/save_tax"] = save_tax

            JAMs_file["integration_info/master_seed"] = master_seed

            if !isnothing(save_functions)
                JAMs_file["integration_info/save_functions"] = save_functions
            end

            if !isnothing(Tplot)

                JAMs_file["integration_info/Tplot"] = Tplot

                JAMs_file["integration_info/fps"] = fps

                JAMs_file["integration_info/plot_functions"] = plot_functions

                JAMs_file["integration_info/plotdim"] = plotdim


            end
            


        end

        #Store similar info in raw_data container
        h5open(joinpath(save_folder_path, raw_data_file_name),"cw") do raw_data_file

            save_raw_metadata!(raw_data_file, system,external_forces,pair_forces, field_forces,local_dofevolvers, global_dofevolvers, field_dofevolvers, integration_tax,dt, t_stop, Tsave,save_tax,master_seed)

        end
   
    end
    if !isnothing(Tsave)
        n_final_save = save_nax[end]
    else
        n_final_save = length(integration_tax)
    end

    if !isnothing(Tplot) 
        cpsO = Observable(current_particle_state)
        cfsO = Observable(current_field_state)
        tO = Observable(0.)
        f, ax = setup_system_plotting(system.sizes,plot_functions, plotdim,cpsO,cfsO,tO,fps,res=res, sbs=sbs)
        

    end
    video_stream = ( !isnothing(record_folder_path) && !isnothing(Tplot) ) ? VideoStream(f, format = format, framerate = fps, visible=true,compression=crf) : nothing
    #Open the files to update over simulation run time
    raw_data_file = !isnothing(Tsave) ? h5open(joinpath(save_folder_path, raw_data_file_name),"r+") : nothing

    #Store frame data here
    frame_group = !isnothing(Tsave) ? create_group(raw_data_file, "frames") : nothing

    #JAMs_file =  !isnothing(Tsave) ? jldopen(joinpath(save_folder_path, JAMs_file_name),"a+") : nothing

    #Loop over time
    #Variable to keep track of the number of frames saved
    frame_counter = 1
    
    try # Catch mechanism to close raw data file in case of an interruption
        @showprogress dt = 1 desc="JAMming in progress..." showspeed=true for (n, t) in pairs(integration_tax)
            
            current_particle_state = threaded_particle_step!(current_particle_state,Next,external_forces, Npair,pair_forces,t, dt, system,cells,rngs_particles)

            if Nfield>0
                for i in eachindex(current_particle_state)
                    p_i = PRow(current_particle_state, i)
                    Forces.contribute_field_forces!(p_i, current_field_state,field_forces, t, dt,system,rngs_particles)
                end
            end

            #Field update loop
            for i in eachindex(current_field_state)

                field_i = current_field_state[i]

                field_i = field_step!(i, field_i,current_field_state,Nfieldu,t,dt, system, rngs_fields)

                current_field_state[i] = field_i
            end

            if !isnothing(Tsave) && !isnothing(save_functions)
                if (n-1)%Tsave==0

                    #Add a new subgroup for this specific frame
                    current_frame_group=create_group(frame_group, string(frame_counter))
    
                    for save_function in save_functions

                        raw_data_file=save_function(current_frame_group,current_particle_state,current_field_state, n, Tsave, t,frame_counter)
                        
                    end

                    frame_counter+=1
                end
            end


            
            #Save the states before the final dof step
            if n==n_final_save
                final_particle_state = deepcopy(current_particle_state)
                final_field_state = deepcopy(current_field_state)



                if !isnothing(Tsave)
                    jldopen(joinpath(save_folder_path, JAMs_final_state_file_name),"a+",iotype=IOStream ) do JAMs_file

                        JAMs_file["SIM"]=SIM(final_particle_state, final_field_state, deepcopy(dt), deepcopy(t_stop), deepcopy(system),true);
                    end
                end

                
            end

            #Only now evolve dofs of every particle
            #Local
            current_particle_state = threaded_dofevolver_step!(current_particle_state,local_dofevolvers,t, dt, system)


            #Or global. The order will first be

            for dofevolver in global_dofevolvers
                current_particle_state,current_field_state = DOFevolvers.evolve_globally!(current_particle_state, current_field_state, system, cells, dt, dofevolver)
            end

            #Perform checks

            current_particle_state = threaded_periodic_bc!(current_particle_state,system)

            cells = update_cells!(cells, current_particle_state)

            #DOF evolver fields
            for i in eachindex(current_field_state)

                field_i = current_field_state[i]

                for dofevolver in field_dofevolvers
                    field_i=evolve_field!(field_i, t, dt, dofevolver)
                end
            end



            
            if !isnothing(Tplot)
                if (n-1)%Tplot==0
                    if !isopen(f.scene)
                        GLMakie_window_closeall()
                        error("Closing the program, because live plotting window is closed.")
                    end
                    notify(cpsO)#[] = current_particle_state
                    notify(cfsO)#[]= current_field_state
                    tO[] = t
                    if !isnothing(video_stream)
                        recordframe!(video_stream)
                    end
                end

            end

        end

        return SIM(deepcopy(current_particle_state), deepcopy(current_field_state), deepcopy(dt), deepcopy(t_stop), deepcopy(system),true);

    catch e

        println("Safely aborting")
        rethrow(e)

    finally 
        if !isnothing(Tsave)
            flush(raw_data_file)
            close(raw_data_file)
        end
        if !isnothing(video_stream)
            save( joinpath(mkpath(record_folder_path),"movie."*format), video_stream)
        end
    end
end

function threaded_particle_step!(current_particle_state,Next,external_forces, Npair,pair_forces,t, dt, system,cells,rngs_particles)
    Threads.@threads for i in eachindex(current_particle_state)

            p_i = PRow(current_particle_state, i)
            init_unwrap!(p_i, t)
            init_f_T!(p_i, t)
            particle_step!(i,p_i,current_particle_state,Next,external_forces, Npair,pair_forces,t, dt, system,cells,rngs_particles)
        end
    return current_particle_state
end

function threaded_dofevolver_step!(current_particle_state,local_dofevolvers,t, dt, system)
    Threads.@threads for i in eachindex(current_particle_state)
        p_i = PRow(current_particle_state, i)
        local_dofevolver_iterate!(p_i, t, dt, local_dofevolvers)
    end
    return current_particle_state
end
function threaded_periodic_bc!(current_particle_state,system)
    Threads.@threads for i in eachindex(current_particle_state)
        p_i = PRow(current_particle_state, i)

        #apply periodic boundary conditions
        if system.Periodic
            periodic!(p_i, system.sizes)
        end
        check_outside_system(i, p_i, system.sizes)
    end
    return current_particle_state
end

#Using foreach is more idiomatic to Julia. 
#However, we want to guarentee that the tuple of forces are unrolled for every length of the tuple.
#Therefore we 'manually' unroll the tuple by using the @generated macro
# function local_dofevolver_iterate!(p_i, t, dt, local_dofevolvers)

#     foreach(dofevolver->evolve_locally!(p_i, t, dt, dofevolver), local_dofevolvers)

#     return p_i
# end
# Use @inline to reduce function call overheads
@generated function local_dofevolver_iterate!(p_i, t, dt, local_dofevolvers::NTuple{N, Any}) where N
    quote
        @nexprs $N k -> @inline DOFevolvers.evolve_locally!(p_i, t, dt, local_dofevolvers[k])
        return p_i
    end
end

# function external_force_iterate!(p_i, t, dt,rngs_particles, system, external_forces)

#     foreach(force->contribute_external_force!(p_i, t, dt,rngs_particles, system, force), external_forces)

#     return p_i
# end

@generated function external_force_iterate!(p_i, t, dt,rngs_particles, system, external_forces::NTuple{N, Any}) where N
    quote
        @nexprs $N k -> @inline Forces.contribute_external_force!(p_i, t, dt,rngs_particles, system, external_forces[k])
        return p_i
    end
    
end

# function pair_force_iterate!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, pair_forces)

#     foreach(force->contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, force),pair_forces)
    
#     return p_i
# end

@generated function pair_force_iterate!(p_i, p_j, dx, dxn, t, dt, rngs_particles, system, pair_forces::NTuple{N, Any}) where N
    quote
        @inline @nexprs $N k -> @inline Forces.contribute_pair_force!(p_i, p_j, dx, dxn, t, dt, rngs_particles, system, pair_forces[k])
        return p_i
    end
end

# function field_force_iterate!(p_i, field_j, field_indices, t, dt,rngs_particles, system, field_forces)
#     foreach(force->contribute_field_force!(p_i, field_j, field_indices, t, dt,rngs_particles, system, force),field_forces)
#     return p_i, field_j
# end

@generated function field_force_iterate!(p_i, field_j, field_indices, t, dt,rngs_particles, system, field_forces::NTuple{N, Any}) where N
    quote
        @nexprs $N k -> @inline Forces.contribute_field_force!(p_i, field_j, field_indices, t, dt,rngs_particles, system, field_forces[k])
        return p_i
    end
end



function particle_step!(i, p_i,current_particle_state, Next,external_forces,Npair,pair_forces,t, dt, system,cells,rngs_particles)
    if Npair>0
        contribute_pair_forces!(i,p_i, current_particle_state ,pair_forces,t, dt, system,cells,rngs_particles)
    end
    if Next>0
        external_force_iterate!(p_i, t, dt,rngs_particles, system, external_forces)
    end
    return p_i
end

function check_outside_system(i,p_i, system_sizes)
    for j in eachindex(p_i.x)
        if p_i.x[j]>system_sizes[j]/2 || p_i.x[j]<-system_sizes[j]/2
            error("JAMs: Particle i = $i is outside simulation box dimension $j. This invalidates cell lists. Suggested fix: make sure particles always stay inside system sizes by increasing the system size of the relevant dimension.")
        end
    end

end

#currently not using current_field_state as there are currently no pair interactions betwen fields, but allow for the implementation.
function field_step!(i, field_i,current_field_state,Nfieldu,t,dt, system, rngs_fields)
    if Nfieldu>0
        for field_updater in system.field_updaters
            field_i=FieldUpdaters.contribute_field_update!(field_i, t, dt, field_updater, rngs_fields)
        end
    end

    return field_i

end


function contribute_field_forces!(p_i, current_field_state,field_forces, t, dt,system, rngs_particles)

    
    for (j, field_j) in pairs(current_field_state)
        field_indices = minimal_image_closest_field_center(p_i.x, field_j.bin_centers, field_j.lbin)

        field_force_iterate!(p_i, field_j, field_indices, t, dt,rngs_particles, system, field_forces)
        

    end
    return p_i, current_field_state
end




function contribute_pair_forces!(i,p_i, current_particle_state, pair_forces, t, dt,system,cells,rngs_particles)


    cell_id_of_p_i = cells.particle_id_to_cell_id[i]

    @inbounds for k in 1:27
        candidate_cell_id = cells.neighbors[k, cell_id_of_p_i]
        candidate_cell_id == 0 && break #if zero: there will be no neigbouring cells anymore so we can stop completely 

        for ind in cells.group_edges[candidate_cell_id]:(cells.group_edges[candidate_cell_id+1]-1)

            j = cells.grouped_particle_ids[ind]

            j == i && continue #If the particle in the candidate cell is particle i, skip the next part

            p_j =  PRow(current_particle_state, j)

            dx = minimal_image_difference(p_i.x, p_j.x, system.sizes, system.Periodic)

            dxn2 = sum(abs2,dx)
            
            if dxn2<=system.rcut_pair_global^2
                dxn = sqrt(dxn2)

                pair_force_iterate!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, pair_forces)

            end
        end
    end
    return p_i
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
                            if !(cell_candidate_index in cell_neighbour_list[:,cell_index])
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
    nbins = ntuple(d -> max(floor(Int, system.sizes[d] / system.rcut_pair_global), 1), 3)
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