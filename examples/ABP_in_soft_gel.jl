
using JAMS

using GLMakie

function relaxation()
    forces = (Forces.repulsive_soft_disk([1,2],[1. 1. ; 1. 1.]),)
    
    dofevolvers = (DOFevolvers.overdamped_xvf([1,2]),DOFevolvers.overdamped_pqT_xyc([1,2]))


    Nconf=1000
    Nact = 100
    N = Nconf+Nact
    ϕ = 0.8
    poly=15e-2
    Rs = vcat( rand(Uniform(1-poly, 1+poly),Nconf),rand(Uniform(1-poly, 1+poly),Nact) )
    display(size(Rs))

    L =  sqrt(pi *sum(Rs.^2) / ϕ)

    types = vcat(ones(Int64,Nconf), ones(Int64,Nact) .*2)

    initial_state = ParticleState([Particles.Polar(id=i,type=types[i], R = Rs[i], x=[rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0],p=normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0])) for i=1:N ]);

    display(L)
    sizes = (L,L,4.);
    print(sizes)

    system =  System(sizes=sizes, initial_particle_state = initial_state,forces = forces, dofevolvers = dofevolvers, Periodic=true,rcut_pair_global=2.5*(1+poly));
    sim = Euler_integrator(system,0.01, 60, Tplot=10,fps=60,plot_functions=(LPlot.disks_type!, LPlot.directors!, LPlot.velocity_vectors!), plotdim=2); 
    return sim;

end
rx = relaxation()
function morse_relaxation(relaxation)

    forces = (Forces.repulsive_soft_disk([1,2],[1. 1. ; 1. 1.]),Forces.morse(1,0.05,2.5))  
    dofevolvers = (DOFevolvers.overdamped_xvf([1,2]),DOFevolvers.overdamped_pqT_xyc([1,2]))


    initial_state = relaxation.final_particle_state
    sizes = relaxation.system.sizes

    
    system = System(sizes=sizes, initial_particle_state = initial_state,forces = forces, dofevolvers = dofevolvers, Periodic=true,rcut_pair_global=10.);

    #Run integration
    #Use plot_disks! for nice visualss
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system,0.05, 5e2, Tplot=10,fps=60,plot_functions=(LPlot.disks_type!,LPlot.velocity_vectors!),plotdim=2)#, record_folder_path=joinpath(homedir(),"soft_gel_19_03-2026_v2"), res= (1000,1000))# plot_velocity_vectors!), plotdim=2 )#, record_folder_path=joinpath(homedir(),"soft_gel_05-01-2026"), res= (1000,1000)); 
    return sim;

end


rx2 = morse_relaxation(rx)  

function simulation(relaxation)

    
    forces = (Forces.self_propulsion(2,0.5), Forces.planar_rotational_noise(ontypes=2, Dr=0.01),Forces.repulsive_soft_disk([1,2],[1. 1. ; 1. 1.]),Forces.morse(1,0.05,2.5))  
    dofevolvers = (DOFevolvers.overdamped_xvf([1,2]),DOFevolvers.overdamped_pqT_xyc([1,2]))

    initial_state = relaxation.final_particle_state
    sizes = relaxation.system.sizes

    system = System(sizes=sizes, initial_particle_state = initial_state,forces = forces, dofevolvers = dofevolvers, Periodic=true,rcut_pair_global=10.);

    sim = Euler_integrator(system,0.05, 5e2, Tplot=10,fps=60,plot_functions=(LPlot.disks_type!,LPlot.velocity_vectors!),plotdim=2)
    return sim;

end

simulation(rx2)
 