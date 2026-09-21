@inline function periodic!(i,p_i,current_particle_state, systemsizes)

    current_particle_state.x[i] = mod.(current_particle_state.x[i] .+ systemsizes ./ 2, systemsizes) .- systemsizes ./ 2
    return current_particle_state
end

@inline function periodic!(i,p_i::Particles.RigidBody,current_particle_state, systemsizes)

current_particle_state.x[i] = mod.(current_particle_state.x[i] .+ systemsizes ./ 2, systemsizes) .- systemsizes ./ 2
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
    return current_particle_state
end

#Initialize unwrapped coordinates to save the user the hassle to set equal to the initial wrapped coordinates
function init_unwrap!(i,current_particle_state, t)

    if t==0
        current_particle_state.xuw[i]= copy(current_particle_state.x[i])
    end

    return current_particle_state
end
#Initialize forces and torques to zero, for easy chaining of sims
function init_f_T!(i,current_particle_state, t)

    if t==0
        current_particle_state.f[i]*=0.
        current_particle_state.T[i]*=0.
    end
    return current_particle_state
end


function minimal_image_closest_bin_center!(field_indices,x, bin_centers,system_sizes,system_Periodic)

    for (i, xi) in pairs(x)

        current_min_ind = 1
        current_min_dis = Inf

        for (j, bcij) in pairs(bin_centers[i])
            dis = abs( bcij - x[i])

            if dis<current_min_dis
                current_min_ind = j
                current_min_dis = dis
            end

        end
        field_indices[i] = current_min_ind

   
    end
    return field_indices

end
#Optimize for field indices
function minimal_image_closest_field_center(x, bin_centers, lbin)

    x_ind = length(bin_centers[1])>1 ? clamp(round(Int, (x[1] - bin_centers[1][2]) / lbin) + 2, 2, length(bin_centers[1])-1) : 1
    y_ind = length(bin_centers[2])>1 ? clamp(round(Int, (x[2] - bin_centers[2][2]) / lbin) + 2, 2, length(bin_centers[2])-1) : 1

    z_ind = length(bin_centers[3])>1 ? clamp(round(Int, (x[3] - bin_centers[3][2]) / lbin) + 2, 2, length(bin_centers[3])-1) : 1

    return  SVector{3, Int64}(x_ind, y_ind, z_ind)

end
#dx is assumed to be pre-allocated
function minimal_image_difference_deprecated!(dx,xi, xj, system_sizes, system_Periodic)
    
    @inbounds for n in eachindex(xi)
        dx[n]=xj[n]-xi[n]
        if system_Periodic
            if dx[n]>system_sizes[n]/2
                dx[n]-=system_sizes[n]
            end
            if dx[n]<=-system_sizes[n]/2
                dx[n]+=system_sizes[n]
            end
        end 
    end
    return dx
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
end
function save_raw_force_data!(file, preamble, force)

    force_name = string(nameof(typeof(force)))

    field_names = fieldnames(typeof(force))

    for field_name in field_names

        name = string(field_name)
        val = getfield(force, field_name)


        file[preamble*force_name*"/"*name] = val

    end

    return file
end

function save_raw_obj_data!(file_group, obj)

    obj_name = string(nameof(typeof(obj)))

    field_names = fieldnames(typeof(obj))

    obj_group = create_group(file_group, obj_name)

    for field_name in field_names

        name = string(field_name)
        val = getfield(obj, field_name)

        if typeof(val)==Bool
            obj_group[name] =string( val)
        else
            obj_group[name] = val
        end

    end

    return file_group
end


function save_raw_metadata!(file, system, integration_tax,dt,t_stop,Tsave,save_tax, master_seed)

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

            save_raw_metadata!(raw_data_file, system, integration_tax,dt, t_stop, Tsave,save_tax,master_seed)

        end
   
    end

    if !isnothing(Tsave)
        n_final_save = save_nax[end]
    else
        n_final_save = length(integration_tax)
    end

    if system.Periodic==false
        @warn ("JAMs: System is set to non-periodic: you should make sure particles always stay in system sizes for correctly working cell lists. The program will catch this by throwing an error if a particle is detected outside the box.")
    end

    system, cells,cell_bin_centers,stencils, lbins = construct_cell_lists!(system)

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

    
    final_particle_state = deepcopy(system.initial_particle_state)
    final_field_state = deepcopy(system.initial_field_state)

    rngs_fields = [Xoshiro(master_seed+i) for i in eachindex(current_field_state)]

    #Assuming number of fields stay constant
    rngs_particles = [Xoshiro(length(current_field_state)+master_seed+i) for i in eachindex(current_particle_state)]

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
            
            current_particle_state = threaded_particle_step!(current_particle_state,Next,external_forces, Npair,pair_forces,t, dt, system,cells,cell_bin_centers,stencils,rngs_particles)

            if Nfield>0
                for i in eachindex(current_particle_state)
                    p_i = current_particle_state[i]
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

            final_particle_state = deepcopy(current_particle_state)
            final_field_state = deepcopy(current_field_state)
            
            #Save the states before the final dof step
            if n==n_final_save



                if !isnothing(Tsave)
                    jldopen(joinpath(save_folder_path, JAMs_final_state_file_name),"a+",iotype=IOStream ) do JAMs_file

                        JAMs_file["SIM"]=SIM(final_particle_state, final_field_state, deepcopy(dt), deepcopy(t_stop), deepcopy(system));
                    end
                end

                
            end

            #Only now evolve dofs of every particle
            #Local
            current_particle_state = threaded_dofevolver_step!(current_particle_state,local_dofevolvers,t, dt, system)


            #Or global. The order will first be

            for dofevolver in global_dofevolvers
                current_particle_state,current_field_state = DOFevolvers.evolve_globally!(current_particle_state, current_field_state, system, cells, stencils, dt, dofevolver)
            end

            #Perform checks

            current_particle_state = threaded_periodic_bc!(current_particle_state,system)

            current_particle_state, cells = update_cells!(current_particle_state, cells, cell_bin_centers,system,lbins)

            #DOF evolver fields
            for i in eachindex(current_field_state)

                field_i = current_field_state[i]

                for dofevolver in field_dofevolvers
                    field_i=evolve_field!(field_i, t, dt, dofevolver)
                end
            end


            cells = update_ghost_cells!(cells,system)

            
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

        return SIM(deepcopy(final_particle_state), deepcopy(final_field_state), deepcopy(dt), deepcopy(t_stop), deepcopy(system));

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

function threaded_particle_step!(current_particle_state,Next,external_forces, Npair,pair_forces,t, dt, system,cells,cell_bin_centers,stencils,rngs_particles)
    Threads.@threads for i in eachindex(current_particle_state)
            p_i = current_particle_state[i]
            init_unwrap!(i,current_particle_state, t)
            init_f_T!(i,current_particle_state, t)
            particle_step!(i, current_particle_state,Next,external_forces, Npair,pair_forces,t, dt, system,cells,cell_bin_centers,stencils,rngs_particles)
        end
    return current_particle_state
end

function threaded_dofevolver_step!(current_particle_state,local_dofevolvers,t, dt, system)
    Threads.@threads for i in eachindex(current_particle_state)

        local_dofevolver_iterate!(i,current_particle_state, t, dt, local_dofevolvers)
    end
    return current_particle_state
end
function threaded_periodic_bc!(current_particle_state,system)
    Threads.@threads for i in eachindex(current_particle_state)
        p_i = current_particle_state[i]

        #apply periodic boundary conditions
        if system.Periodic
            periodic!(i,p_i,current_particle_state, system.sizes)
        end
        check_outside_system(i,current_particle_state, system.sizes)
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
@generated function local_dofevolver_iterate!(i,current_particle_state, t, dt, local_dofevolvers::NTuple{N, Any}) where N
    quote
        @nexprs $N k -> @inline DOFevolvers.evolve_locally!(i,current_particle_state, t, dt, local_dofevolvers[k])
        return current_particle_state
    end
end

# function external_force_iterate!(p_i, t, dt,rngs_particles, system, external_forces)

#     foreach(force->contribute_external_force!(p_i, t, dt,rngs_particles, system, force), external_forces)

#     return p_i
# end

@generated function external_force_iterate!(i,current_particle_state, t, dt,rngs_particles, system, external_forces::NTuple{N, Any}) where N
    quote
        @nexprs $N k -> @inline Forces.contribute_external_force!(i,current_particle_state, t, dt,rngs_particles, system, external_forces[k])
        return current_particle_state
    end
    
end

# function pair_force_iterate!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, pair_forces)

#     foreach(force->contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, force),pair_forces)
    
#     return p_i
# end

@generated function pair_force_iterate!(i,p_i, p_j, current_particle_state, dx, dxn, t, dt, rngs_particles, system, pair_forces::NTuple{N, Any}) where N
    quote
        @inline @nexprs $N k -> @inline Forces.contribute_pair_force!(i,p_i, p_j, current_particle_state, dx, dxn, t, dt, rngs_particles, system, pair_forces[k])
        return current_particle_state
    end
end

# function field_force_iterate!(p_i, field_j, field_indices, t, dt,rngs_particles, system, field_forces)
#     foreach(force->contribute_field_force!(p_i, field_j, field_indices, t, dt,rngs_particles, system, force),field_forces)
#     return p_i, field_j
# end

@generated function field_force_iterate!(p_i, field_j, field_indices, t, dt,rngs_particles, system, field_forces::NTuple{N, Any}) where N
    quote
        @nexprs $N k -> @inline Forces.contribute_field_force!(p_i, field_j, field_indices, t, dt,rngs_particles, system, field_forces[k])
        return current_particle_state
    end
end



function particle_step!(i,current_particle_state,Next,external_forces,Npair,pair_forces,t, dt, system,cells,cell_bin_centers,stencils,rngs_particles)
    if Npair>0
        contribute_pair_forces!(i, current_particle_state,pair_forces,t, dt, system,cells,stencils,rngs_particles)
    end
    if Next>0
        external_force_iterate!(i, current_particle_state, t, dt,rngs_particles, system, external_forces)
    end
    return current_particle_state
end

function check_outside_system(i, current_particle_state, system_sizes)
    p_i = current_particle_state[i]
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




function contribute_pair_forces!(i, current_particle_state, pair_forces, t, dt,system,cells,stencils,rngs_particles)
    @inbounds p_i = current_particle_state[i]
    for stencil in stencils
        candidate_cell_inds = p_i.ci .+ stencil

        @inbounds for n in cells[candidate_cell_inds...]
            if i!=n
                p_j = current_particle_state[n]

                dx = minimal_image_difference(p_i.x, p_j.x, system.sizes, system.Periodic)

                dxn2 = sum(abs2,dx)
                
                if dxn2<=system.rcut_pair_global^2
                    dxn = sqrt(dxn2)

                    pair_force_iterate!(i,p_i, p_j,current_particle_state, dx, dxn, t, dt,rngs_particles, system, pair_forces)

                end

            end
        end
    end
end



function construct_cell_list_centers(L,rcut_pair_global)

    #Choose bin size slight larger if lbin does not divide box length
    nbin = floor(Int64,L/rcut_pair_global)

    #In case user sets irrelevant dimension smaller than rcut_pair_global

    if nbin<1
        nbin=2 #very important, must have at least 2 real bins besides the ghost cells 
    end
    lbin = L/nbin


    bin_centers = Float64[-L-rcut_pair_global]
    #0.1 factor to make sure the end point is inluded
    bin_centers = append!(bin_centers,range(start=-(L-lbin)/2, stop=(L-lbin)/2+0.1*lbin, step=lbin))
    bin_centers = append!(bin_centers,L+rcut_pair_global)

    return bin_centers, lbin

end

function construct_cell_lists!(system)

    dim = length(system.sizes)

    Lx = system.sizes[1]
    Ly = system.sizes[2]
    Lz = system.sizes[3]

    #First and last bin center should be far outside the simulation box so a particle is never associated in a ghost cell
    #via the minimal_image_closest_bin_center function. Note that system spans from -Li/2 to Li/2 so -Li should be
    #safely outside reach both when using periodic boundary conditions
    #or when using a finite system, because particle outside box will raise and error and  abort program
    #In case system size is entered as Int64, without FLoat64 in front, the array will be typed as Int64,
    # not allowing appending the  float values in the next line. Therefore declare type as Float64[].
    # +lbin incase lbin is larger than one of the system sizes
    x_bin_centers,x_lbin = construct_cell_list_centers(Lx,system.rcut_pair_global)
    y_bin_centers, y_lbin = construct_cell_list_centers(Ly,system.rcut_pair_global)
    z_bin_centers,z_lbin = construct_cell_list_centers(Lz,system.rcut_pair_global)
    nx = length(x_bin_centers)
    ny = length(y_bin_centers)
    nz = length(z_bin_centers)

    cell_bin_centers = [x_bin_centers, y_bin_centers, z_bin_centers]

    cells = reshape([Int64[] for i in 1:nx*ny*nz],nx,ny,nz)
    for (i,p_i) in pairs(system.initial_particle_state)

        cell_indices=@MVector zeros(Int,length(p_i.x))
        cell_indices = minimal_image_closest_bin_center!(cell_indices,p_i.x, cell_bin_centers,system.sizes,system.Periodic)

        cells[cell_indices...]= append!(cells[cell_indices...],p_i.id)
        system.initial_particle_state.ci[i] = cell_indices

    end
    stencils = [ @SVector [ni, nj, nk] for ni in -1:1 for nj in -1:1 for nk in -1:1]

    lbins=(x_lbin, y_lbin, z_lbin)

    cells = update_ghost_cells!(cells,system)

    return system, cells, cell_bin_centers, stencils, lbins
end

@inbounds function find_new_bin_location(p_i, cell_bin_centers,system,lbins)

    moved=false
    #Collect old locations to preallocate for new one
    x_ind = length(cell_bin_centers[1])>1 ? clamp(round(Int, (p_i.x[1] - cell_bin_centers[1][2]) / lbins[1]) + 2, 2, length(cell_bin_centers[1])-1) : 1
    moved = moved || (x_ind != p_i.ci[1])

    y_ind = length(cell_bin_centers[2])>1 ? clamp(round(Int, (p_i.x[2] - cell_bin_centers[2][2]) / lbins[2]) + 2, 2, length(cell_bin_centers[2])-1) : 1
    moved = moved || (y_ind != p_i.ci[2])

    z_ind = length(cell_bin_centers[3])>1 ? clamp(round(Int, (p_i.x[3] - cell_bin_centers[3][2]) / lbins[3]) + 2, 2, length(cell_bin_centers[3])-1) : 1
    moved = moved || (z_ind != p_i.ci[3])

    return SVector{3, Int64}(x_ind, y_ind, z_ind), moved

end



function update_cells!(current_particle_state, cells, cell_bin_centers,system,lbins)

    #Updating the cell list is not threadsafe, hence put it outside the threaded loop
    #The cell list update must come *after* the DOF evolving of particles!

    #new_bin_location = @MVector  zeros(Int64, length(current_particle_state[1].ci))
    for i in eachindex(current_particle_state)

        p_i = current_particle_state[i]
        #initialize with the old bin location
        #copyto!(new_bin_location, p_i.ci)
        new_bin_location,moved=find_new_bin_location(p_i, cell_bin_centers,system,lbins)
        if moved

            cell_view = @views cells[p_i.ci...]
            ind = findfirst(x -> x == p_i.id, cell_view)
            # #Swap and pop

            cells[p_i.ci...][end], cells[p_i.ci...][ind] = cells[p_i.ci...][ind], cells[p_i.ci...][end]
            pop!(cells[p_i.ci...])
            
            #add to correct lists
            current_particle_state.ci[i]= new_bin_location

            push!(cells[new_bin_location...],p_i.id)
        end
    end
    return current_particle_state, cells
end



function update_ghost_cells!(cells,system)
    if system.Periodic
        
        dims = length(system.sizes)
        if dims==2

            cells[1,:] .= @view cells[end-1,:]
            cells[end,:] .= @view cells[2,:]

            cells[:,1] .= @view cells[:,end-1]
            cells[:,end] .= @view cells[:,2]
        end

        if dims==3

            #The bottom works, but one can do it much more easily: imagine e.g. 5 by 5 cube and see how each points gets correctly
            # assigned by simply doing the following steps IN ORDER

            cells[1, :, :]   .= @view cells[end-1, :, :]
            cells[end, :, :] .= @view cells[2, :, :]

            cells[:, 1, :]   .= @view cells[:, end-1, :]
            cells[:, end, :] .= @view cells[:, 2, :]

            cells[:, :, 1]   .= @view cells[:, :, end-1]
            cells[:, :, end] .= @view cells[:, :, 2]

    
        end
    end
    return cells
end
