
using JAMS

using GLMakie

function simulation()

    forces = (Forces.self_align_with_v(1,0.,false),Forces.self_propulsion(1,0.2),Forces.rotational_noise_3d(ontypes=1,Dr=0.01),Forces.repulsive_soft_disk(1,1.),)

    dofevolvers = (DOFevolvers.overdamped_xvf(1),DOFevolvers.overdamped_pqT(1))

    N=100
    ϕ = 0.1
    poly=15e-2
    Rs =rand(Uniform(1-poly, 1+poly),N)
    display(size(Rs))

    L =  ( 4/3 * pi *sum(Rs.^3) / ϕ)^(1/3)

    initial_state = ParticleState([ Particles.Polar(id=i,type=1,R=Rs[i],  x=rand(Uniform(-L/2, L/2),3) , p=normalize(rand(Normal(0, 1),3))) for i=1:N ])


    display(L)
    sizes = (L,L,L);
    print(sizes)


    system = System(sizes=sizes, initial_particle_state = initial_state,forces = forces, dofevolvers = dofevolvers, Periodic=true,rcut_pair_global=2.5*(1+poly));
    
    sim = Euler_integrator(system,0.01, 1e4,Tplot=10,fps=60,plot_functions=(LPlot.sized_points!,LPlot.directors!,LPlot.velocity_vectors!),plotdim=3); 
    return sim;

end


sim = simulation()



