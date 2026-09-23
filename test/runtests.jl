using JAMS
using Test

@testset "JAMS.jl" begin
    
    @testset "Micro headless simulation run" begin

    forces = (Forces.self_align_with_v(1,0.1,false),Forces.self_propulsion(1,0.2),Forces.planar_rotational_noise(ontypes=1,Dr=.01),Forces.repulsive_soft_disk(1,1.),)
    dofevolvers = (DOFevolvers.overdamped_xvf(1),DOFevolvers.overdamped_pqT_xyc(1))
    N=10
    ϕ = 1.0
    poly=15e-2
    Rs =rand(Uniform(1-poly, 1+poly),N)
    L =  sqrt(pi *sum(Rs.^2) / ϕ)
    initial_state = ParticleState([ Particles.Polar(id=i,type=1,R=Rs[i], x=[rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0],p=normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0])) for i=1:N ])

    sizes = (L,L,2.);

    system = System(sizes=sizes, initial_particle_state = initial_state,forces = forces, dofevolvers = dofevolvers, Periodic=true,rcut_pair_global=2.5*(1+poly));
    
    sim = Euler_integrator(system,0.01, 1); 
        


    end



end
