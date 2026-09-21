using JAMS

using GLMakie

function relaxation()
    forces = (Forces.repulsive_soft_disk([1,2],[1 2 ; 2 2]),)

    dofevolvers = (DOFevolvers.overdamped_xvf(1),DOFevolvers.overdamped_pqT_xyc(1))

    aspect = 1
    Lx = 50. *aspect
    Ly = Lx /aspect^2
    phi = 1.3
    poly=15e-9
    l = 1.5
    xs, ys = Initial.box(l, Lx, Ly)
    Nb = length(xs)
    Rb =rand(Uniform(1-poly, 1+poly),Nb)

    A = Lx*Ly - pi *sum(Rb.^2)/2
    Nint = trunc(Int64,cld(A*phi,pi))
    Rint = rand(Uniform(1-poly, 1+poly),Nint)
    display(Nint+Nb)

    pd = [ Particles.Polar(id=i,type=1,R=Rint[i], x=[rand(Uniform(-Lx/2 +4, Lx/2 - 4)) , rand(Uniform(-Ly/2 +4, Ly/2 - 4)),0],p=normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0])) for i=1:Nint ]

    initial_state = ParticleState( append!(pd,[ Particles.Polar(id=i+Nint,type=2,R=Rb[i], x=[xs[i] , ys[i],0],p=normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0])) for i=1:Nb ]))
    

    sizes = (Lx+4, Ly+4,1.);
    print(sizes)

    system = System(sizes=sizes, initial_particle_state = initial_state,forces = forces, dofevolvers = dofevolvers, Periodic=true,rcut_pair_global=2.5*(1+poly));

    sim = Euler_integrator(system,0.05, 1e3,Tplot=100,fps=60,plot_functions=(LPlot.disks_v_orientation!,LPlot.directors!),plotdim=2); 
    return sim;

end
rx=relaxation()


function sa_step(rx)
    forces = (Forces.self_align_with_v(1,0.1,true),Forces.self_propulsion(1,0.01),Forces.planar_rotational_noise(ontypes=1,Dr=0.01),Forces.repulsive_soft_disk([1,2],[1 2 ; 2 2]),)

    dofevolvers = (DOFevolvers.overdamped_xvf(1),DOFevolvers.overdamped_pqT_xyc(1))


    
    initial_state = deepcopy(rx.final_particle_state)
    display(length(initial_state))


    system = System(sizes=rx.system.sizes, initial_particle_state = initial_state,forces = forces, dofevolvers = dofevolvers, Periodic=false,rcut_pair_global=rx.system.rcut_pair_global);

    sim = Euler_integrator(system,0.05, 1e4,Tplot=10,fps=60,plot_functions=(LPlot.disks_v_orientation!,LPlot.directors!),plotdim=2); 
    return sim;

end

sa=sa_step(rx)
