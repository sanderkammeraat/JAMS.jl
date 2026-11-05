
using JAMS

using GLMakie

function simulation_3d()

    compressed = 0.6
    forces = (Forces.self_propulsion(1,0.1),Forces.translational_noise(1,1e-3*[1,1,1]),Forces.polymer_repulsive_soft_disk([1,2],ones(2,2)),Forces.polymer_harmonic_stretch([1,2],false,ones(2,2), compressed*ones(2,2)))

    dofevolvers = (DOFevolvers.overdamped_xvf([1,2]),
        #DOFevolvers.overdamped_polymer_pairdis([1,2],2*ones(2,2), 1.2*ones(2,2) ,true),
        DOFevolvers.polymer_p_set([1,2]))
    N=5
    ϕ = 0.01
    poly=15e-2
    Rs =rand(Uniform(1-poly, 1+poly),N)
    display(size(Rs))

    L =  ( 2*compressed * 4/3 * pi * N  / ϕ)^(1/3)

    initial_state =Particles.PolarPolymer[]

    for i = 1:N
        id1 = 2i -1 
        id2 = 2i
        x1 = rand(Uniform(-L/2, L/2),3)
        p1 = normalize(rand(Normal(0, 1),3))
        push!(initial_state, Particles.PolarPolymer(id = id1,type=1, id_in_pol = 2, pol_id = id1, pol_N = 2, x=x1 , p=p1 ))
        x2 = x1 .- p1 * compressed
        push!(initial_state, Particles.PolarPolymer(id = id2, type=2, id_in_pol = 1, pol_id = id1, pol_N = 2, x=x2 , p=p1 ))
    end

    initial_state = ParticleState(initial_state)

    display(L)
    sizes = (L,L,L);
    print(sizes)


    system = System(sizes=sizes, initial_particle_state = initial_state,forces = forces, dofevolvers = dofevolvers, Periodic=true,rcut_pair_global=2.5*(1+poly));
    
    sim = Euler_integrator(system,0.01, 1e4,Tplot=10,fps=60,plot_functions=(
        LPlot.polymers_3d!,
        #LPlot.type_sized_points!,
    LPlot.directors!,LPlot.velocity_vectors!),plotdim=3); 
    return sim;

end


sim = simulation_3d()


function simulation_2d()

    compressed = 0.6

    forces = (Forces.self_propulsion(1,0.3),Forces.translational_noise(1,1e-3*[1,1,0.]),Forces.polymer_repulsive_soft_disk([1,2],ones(2,2)),Forces.polymer_harmonic_stretch([1,2],false,3. *ones(2,2), compressed*ones(2,2)))

    dofevolvers = (DOFevolvers.overdamped_xvf([1,2]),
                #DOFevolvers.overdamped_polymer_pairdis([1,2],2*ones(2,2), 1.5*ones(2,2) ,true),
                DOFevolvers.polymer_p_set([1,2]))

    N=1000
    ϕ = 0.6

    L =  ( 2. * pi * N * compressed / ϕ)^(1/2)

    initial_state =Particles.PolarPolymer[]

    for i = 1:N
        id1 = 2i -1 
        id2 = 2i
        x1 = [rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0]
        p1 = normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0])
        push!(initial_state, Particles.PolarPolymer(id = id1,type=1, id_in_pol = 2, pol_id = id1, pol_N = 2, x=x1 , p=p1 ))
        x2 = x1 .- p1 *compressed
        push!(initial_state, Particles.PolarPolymer(id = id2, type=2, id_in_pol = 1, pol_id = id1, pol_N = 2, x=x2 , p=p1 ))
    end

    initial_state = ParticleState(initial_state)

    display(L)
    sizes = (L,L,2.);
    print(sizes)


    system = System(sizes=sizes, initial_particle_state = initial_state,forces = forces, dofevolvers = dofevolvers, Periodic=true,rcut_pair_global=2.5);
    
    sim = Euler_integrator(system,0.01, 1e4,Tplot=10,fps=60,plot_functions=(
        LPlot.polymers!,
        #LPlot.type_sized_points!,
    LPlot.directors!,LPlot.velocity_vectors!),plotdim=2); 
    return sim;

end


sim = simulation_2d()




