
using JAMS

using GLMakie

function simulation()

    forces = (Forces.self_align_with_v(1,0.0,false),Forces.self_propulsion(1,0.3),Forces.planar_rotational_noise(ontypes=1,Dr=.01),Forces.repulsive_soft_disk(1,1.),)

    dofevolvers = (DOFevolvers.overdamped_xvf(1),DOFevolvers.overdamped_pqT_xyc(1))

    N=1000
    ϕ = 1.0
    poly=15e-2
    Rs =rand(Uniform(1-poly, 1+poly),N)
    display(size(Rs))

    L =  sqrt(pi *sum(Rs.^2) / ϕ)

    initial_state = ParticleState([ Particles.Polar(id=i,type=1,R=Rs[i], x=[rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0],p=normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0])) for i=1:N ])


    display(L)
    sizes = (L,L,2.);
    print(sizes)


    system = System(sizes=sizes, initial_particle_state = initial_state,forces = forces, dofevolvers = dofevolvers, Periodic=true,rcut_pair_global=2.5*(1+poly));

    
    sim = Euler_integrator(system,0.01, 1e4,Tplot=10,fps=60,plot_functions=(LPlot.disks_v_orientation!,LPlot.directors!),plotdim=2); 
    return sim;

end


sim = simulation()

