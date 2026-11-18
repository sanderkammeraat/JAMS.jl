using JAMS

using GLMakie

function custom_sized_points!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ p_i.id[1] for p_i in $cpsO])
        
    z = @lift([p_i.x[3] for p_i in $cpsO])

    R = @lift([p_i.R[1] for p_i in $cpsO])

    meshscatter!(ax,x,y,z, color=:white, markersize =R, transparency=true, alpha=0.2)#, colormap=(:viridis,0.2))
    
    return ax
end




function custom_directors!(f,ax, cpsO, cfsO)

    x = @lift([ Point3f( p_i.x[1],p_i.x[2],p_i.x[3]) for p_i in $cpsO])
    p = @lift([ Point3f( p_i.p[1],p_i.p[2],p_i.p[3]) for p_i in $cpsO])
    c = @lift([ angle(p_i.p[1]+1im*p_i.p[2]) for p_i in $cpsO])

    arrows3d!(ax, x, p , color=c,  colormap=:hsv,colorrange=(-pi,pi), lengthscale = 2)

    return ax

end



function simulation()

    kpar = 1
    kper = 0
    a = 2.5
    D = 0.2
    epsilon = 1.

    forces = (Forces.pair_polar_alignment(1, 2., .3),Forces.morse(1,0.08,2.5),
    Forces.self_propulsion(1,0.02),Forces.rotational_noise_3d(ontypes=1,Dr=0.001),Forces.repulsive_soft_disk(2,1.),Forces.pairAN(1,true,false,2., kpar, kper, 0.1))

    dofevolvers = (DOFevolvers.overdamped_xvf(1),DOFevolvers.overdamped_pqT(1))

    N=1000
    ϕ = 0.02
    poly=15e-2
    Rs =rand(Uniform(1-poly, 1+poly),N)
    display(size(Rs))

    L =  ( 4/3 * pi *sum(Rs.^3) / ϕ)^(1/3)

    initial_state = ParticleState([ Particles.Polar(id=i,type=1,R=Rs[i],  x=rand(Uniform(-L/2, L/2),3) , p=normalize(rand(Normal(0, 1),3))) for i=1:N ])


    display(L)
    sizes = (L,L,L);
    print(sizes)


    system = System(sizes=sizes, initial_particle_state = initial_state,forces = forces, dofevolvers = dofevolvers, Periodic=true,rcut_pair_global=2.5*(1+poly));
    
    sim = Euler_integrator(system,0.05, 5e4,Tplot=10,fps=60,plot_functions=(custom_sized_points!,custom_directors!,LPlot.velocity_vectors!),plotdim=3); 


    return sim;

end

sim = simulation()

