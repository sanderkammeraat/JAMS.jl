using JAMS

using GLMakie

function relaxation()
    forces = (Forces.repulsive_soft_disk([1,2],[1 2 ; 2 2]),)

    dofevolvers = (DOFevolvers.overdamped_xvf(1),DOFevolvers.overdamped_pqT_xyc(1))

    phi = 1.1
    poly=15e-2
    l = 1.5
    Rin = 0.
    Rout = 40.
    xin, yin = Initial.ring(l, Rin)
    xout, yout = Initial.ring(l,Rout)

    xb = vcat(xin, xout)
    yb = vcat(yin, yout)

    Nb = length(xb)
    Rb =rand(Uniform(1-poly, 1+poly),Nb)

    A =  pi*Rout^2 - pi*Rin^2 - pi *sum(Rb.^2)/2
    Nint = trunc(Int64,cld(A*phi,pi))

    Rint = rand(Uniform(1-poly, 1+poly),Nint)
    display(Nint+Nb)

    xint, yint = Initial.random_in_disk(Nint, Rout-2, R_in=Rin)

    pd = [ Particles.Polar(id=i,type=1,R=Rint[i], x=[xint[i] , yint[i],0],p=normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0])) for i=1:Nint ]

    initial_state = ParticleState( append!(pd,[ Particles.Polar(id=i+Nint,type=2,R=Rb[i], x=[xb[i] , yb[i],0],p=normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0])) for i=1:Nb ]))
    
    Lx = 2*Rout + 2
    Ly = 2*Rout + 2
    sizes = (Lx+4, Ly+4,1.);
    print(sizes)

    system = System(sizes=sizes, initial_particle_state = initial_state,forces = forces, dofevolvers = dofevolvers, Periodic=true,rcut_pair_global=2.1*(1+poly));

    sim = Euler_integrator(system,0.05, 1e3,Tplot=100,fps=60,plot_functions=(LPlot.disks_v_orientation!,LPlot.directors!),plotdim=2); 
    return sim;

end
rx=relaxation()


function sa_step(rx)    
    forces = (Forces.self_align_with_v(1,0.1,true),Forces.self_propulsion(1,0.001),Forces.planar_rotational_noise(ontypes=1,Dr=0.001),Forces.repulsive_soft_disk([1,2],[1 1 ; 1 1]),)

    dofevolvers = (DOFevolvers.overdamped_xvf(1),DOFevolvers.overdamped_pqT_xyc(1))


    
    initial_state = deepcopy(rx.final_particle_state)
    display(length(initial_state))


    system = System(sizes=rx.system.sizes, initial_particle_state = initial_state,forces = forces, dofevolvers = dofevolvers, Periodic=false,rcut_pair_global=rx.system.rcut_pair_global);

    sim = Euler_integrator(system,0.05, 1e4,Tplot=10,fps=60,plot_functions=(LPlot.disks_v_orientation!,LPlot.directors!),plotdim=2); 
    return sim;

end

sa=sa_step(rx)
