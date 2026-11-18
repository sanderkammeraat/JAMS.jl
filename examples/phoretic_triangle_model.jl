using JAMS
using GLMakie
using Interpolations
cvals = log10.([ 4.98472E-6,   1.00375E-5,  5.02894E-5, 9.94311E-5, 3.0171E-4, 4.98163E-4, 6.55346E-4, 8.48719E-4, 0.00102, 0.00503 ] ./ 5e-4)
vvals = [ 0.02162,      0.3501,      0.95642,    1.88254,    2.98665,   4.20652,     4.99807,    3.3817,    0.80448, 0.02162 ]

function make_intpd(cvals, vvals)
    itp = linear_interpolation(cvals, vvals; extrapolation_bc = 0.0)
    return x -> itp(log10(x))
end
intpd = make_intpd(cvals, vvals)

#sim units
# [t] = s
# [L] =  um

function simulation(intpd)

    R = 5.
    A = 30
    a = 5*0.4
    dofevolvers =  (DOFevolvers.overdamped_xvf(1), DOFevolvers.overdamped_pqT_xyc(1),
     DOFevolvers.overdamped_CCvCf(1), DOFevolvers.rigid_xsrel_ysrel(1))

    phi = 0.1
    L = 500.
    N = floor(Int64, L^2 * phi / ( π * R^2))

    display(N)
    Lx=L
    Ly=L
    sizes = (Lx,Ly,1.);
    lbin = 2.   
    display(lbin)


    #Consumption is tied to dirac delta, which has inverse bin size^2 dimensions which need to be compensated
    # One R comes from lbin, the other for scaling the magnitudes of the active forces due to the speed, which should not affect the torque
    #xk, yk  = Initial.circular_cap(pi/5, R, lbin)

    J = 2.0
    beta= 120
    kd = 300
    r = 1.5
    xk = [0.]

    x_vs = [ -4/3*sqrt(3)*r ,-1/3*sqrt(3)*r   , 2/3*sqrt(3)*r ,  2/3*sqrt(3)*r , 2/3*sqrt(3)*r , -1/3*sqrt(3)*r]
    y_vs = [ 0, r, 2r, 0, -2r, -r,]

    triangle = Initial.construct_populated_polygon(x_vs, y_vs, 2)

    xsrel0 = triangle.xsrel0
    ysrel0 = triangle.ysrel0
    rs = r .* ones(length(xsrel0)) 
    #initial_state = ParticleState([ Particles.PolarKernel(id=i,type=1,R =R, x = [rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0],p = normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]), kernel_x = kernel_x, kernel_y = kernel_y) for i=1:N ] );
    initial_state = ParticleState([ Particles.PolarShapeKernel(id=i,type=1,R =R, x = [rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0],p = normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]), xsrel0=xsrel0, ysrel0=ysrel0, Rs=rs) for i=1:N ] );
    forces = (Forces.exp_repulsion(1, A , a),
              Forces.shape_repulsive_soft_disk(1,1.),
              Forces.planar_rotational_noise(ontypes=1,Dr=.001),
              Forces.field_propulsion_interpolation(ontypes = 1, consumption = kd/(length(xk)*lbin^2), interpolation_function = intpd),
              Forces.align_with_gradC(1,J),
              Forces.phoresis(1,-1*beta))

    x_bin_centers = [-Lx-lbin]
    x_bin_centers = append!(x_bin_centers,range(start=-Lx/2, stop=Lx/2-lbin/2, step=lbin).+lbin/2)
    x_bin_centers = append!(x_bin_centers,Lx+lbin/2)

    y_bin_centers = [-Ly-lbin]
    y_bin_centers = append!(y_bin_centers,range(start=-Ly/2, stop=Ly/2-lbin/2, step=lbin).+lbin/2)
    y_bin_centers = append!(y_bin_centers,Ly)
    z_bin_centers = [0.]
    bin_centers = [x_bin_centers, y_bin_centers,z_bin_centers]

    C = 1 .*ones(length(x_bin_centers), length(y_bin_centers))

    display(size(C))
    
    initial_field_state=[Fields.Scalar2d(1,1,bin_centers,C, C.*0, C.*0, lbin)]
    Dc = 400
    kr = 0.3
    field_updaters = (FieldUpdaters.Diffusion(1,Dc), FieldUpdaters.Relax(1,1.,kr), FieldUpdaters.GhostSetPeriodic(1))


    system = System(sizes = sizes, initial_particle_state = initial_state,initial_field_state =initial_field_state , forces = forces, field_updaters =field_updaters, dofevolvers=dofevolvers, Periodic = true, rcut_pair_global=10*R);


    sim = Euler_integrator(system,0.002,1e2,  Tplot=10,
     plot_functions=(LPlot.field_magnitude!, LPlot.disks_orientation!, LPlot.shape_disks_orientation!,LPlot.directors!,), plotdim=2)#

end

sim = simulation(intpd)

