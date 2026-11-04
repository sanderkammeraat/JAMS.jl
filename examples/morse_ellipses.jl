using JAMS

using GLMakie


function relax()

    epsilon = 1.
    beta = 2.
    forces = (Forces.repulsive_soft_ellipse(1, epsilon, beta),)

    #For nondeformable ellipses, use this:
    #dofevolvers = (DOFevolvers.overdamped_xvf(1),DOFevolvers.overdamped_pqT_xyc(1))
    #For deformable ellipses use this
    dofevolvers = (DOFevolvers.overdamped_xvf(1), DOFevolvers.overdamped_deform(1,  false) )

    N=300
    R0 = 1.
    aspect_ratio0 = 1.2

    sigma_edge = R0 * sqrt(aspect_ratio0)
    sigma_side = R0 / sqrt(aspect_ratio0)

    lmda_major_0 = sigma_edge
    lmda_minor_0 = sigma_side
    semi_major = lmda_major_0

    semi_minor = lmda_minor_0
    K =1. 
    mu = 0.5
    tau = 10.
    mu_alpha=0.01
    ϕ = 0.8
    area = pi * lmda_major_0 * lmda_minor_0
    L =  sqrt(N*area/ ϕ)

    thetas = rand(Uniform(-pi, pi), N)
    sins = sin.(thetas)
    coss = cos.(thetas)
    Lambdas = [   [ coss[i]^2*lmda_major_0+sins[i]^2*lmda_minor_0  sins[i]*coss[i]*(lmda_major_0-lmda_minor_0) ; sins[i]*coss[i]*(lmda_major_0-lmda_minor_0) sins[i]^2*lmda_major_0+coss[i]^2*lmda_minor_0 ] for i=1:N ]

    initial_state = ParticleState([ Particles.GBEllipse(id=i,type=1,R=R0,zeta_R =4*R0^2/3  , x=[rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0], p=[ cos(thetas[i]), sin(thetas[i]), 0.],Lambda = Lambdas[i], Lambda0 = Lambdas[i], K = K, mu = mu, tau= tau, mu_alpha=mu_alpha,
     aspect_ratio0 = aspect_ratio0 , semi_major = semi_major, semi_minor = semi_minor, semi_major0 = semi_major, semi_minor0 = semi_minor ) for i=1:N ])

    display(L)
    sizes = (L,L,2.);
    print(sizes)

    system = System(sizes=sizes, initial_particle_state = initial_state,forces = forces, dofevolvers = dofevolvers, Periodic=true,rcut_pair_global=6.);
    
    sim = Euler_integrator(system,0.01, 100,Tplot=10,fps=60,plot_functions=(LPlot.ellipses!,LPlot.directors!),plotdim=2); 
    return sim;

end

rx = relax()

function simulation(rx)

    a = 2.5
    D = 0.2
    epsilon = 1.
    forces = (Forces.self_propulsion(ontypes = 1,v0 = 0.1), Forces.planar_rotational_noise(ontypes=1,Dr=.01),Forces.morse_ellipse(1, a, D, epsilon))    #For nondeformable ellipses, use this:
    
    #For nondeformable ellipses, use this:
    #dofevolvers = (DOFevolvers.overdamped_xvf(1),DOFevolvers.overdamped_pqT_xyc(1))
    #For deformable ellipses use this
    dofevolvers = (DOFevolvers.overdamped_xvf(1), DOFevolvers.overdamped_deform(1,  false) )


    initial_state = rx.final_particle_state
    sizes = rx.system.sizes

    
    system = System(sizes=sizes, initial_particle_state = initial_state,forces = forces, dofevolvers = dofevolvers, Periodic=true,rcut_pair_global=6.);

    sim = Euler_integrator(system,0.01, 1e3,Tplot=30,fps=60,plot_functions=(LPlot.ellipses!,LPlot.directors!),plotdim=2); 
    return sim;

end

sim = simulation(rx);
