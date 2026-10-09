
using JAMS

using GLMakie

function simulation()


    kbend= 3.
    kpar = -1
    N_in_pol = 10
    p=0.3
    f_eq_stretch_force = .75
    krep = 1.
    kbind = 3.

    #forces = (Forces.polymer_pairAN(1,false, false, false, 1.5, kpar, 0., p),Forces.polymer_repulsive_soft_disk(1,krep), Forces.polymer_harmonic_bend(1, kbend),  Forces.polymer_harmonic_stretch(1, false, kbind, f_eq_stretch_force),Forces.translational_noise(1,1e-5*[1,1,0]))
    forces = (Forces.polymer_pair_polar_nematic(1,false, 1.2,p),Forces.polymer_repulsive_soft_disk(1,krep), Forces.polymer_harmonic_bend(1, kbend),  Forces.polymer_harmonic_stretch(1, false, kbind, f_eq_stretch_force),Forces.translational_noise(1,1e-5*[1,1,0]))


    dofevolvers = (DOFevolvers.overdamped_xvf(1),DOFevolvers.polymer_p_set(1))

    pf = 0.9
    R = 1
    L_0 = 80

    x, y, radii, pol_ids, ids_in_pol, L, Npols = Initial.stacked_polymers_at_angle(N_in_pol, R, pf, f_eq_stretch_force, L_0, random_polarity = false)

    initial_state = ParticleState([Particles.PolarPolymer(id=id,type=1,pol_id=pol_ids[id],id_in_pol=ids_in_pol[id],pol_N=N_in_pol, R = radii[id], x=[x[id] , y[id],0],p=normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0])) for id=1:Npols*N_in_pol]);

    display(Npols*N_in_pol)

    sizes = (L,L,2.);


    
    system =  System(sizes=sizes, initial_particle_state = initial_state,forces = forces, dofevolvers = dofevolvers, Periodic=true,rcut_pair_global=6.)

    sim = Euler_integrator(system, 0.01, 1e4, fps=60, Tplot=10, plot_functions =  (LPlot.polymers!, LPlot.nematic_directors! ,LPlot.velocity_vectors!));#, record_folder_path = pwd()); 
    return sim;

end

simulation()