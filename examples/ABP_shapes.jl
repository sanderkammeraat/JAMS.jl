using JAMS
using GLMakie


function custom_directors!(f,ax, cpsO, cfsO)


    x = @lift([ Point3f( p_i.x[1],p_i.x[2],p_i.x[3]) for p_i in $cpsO])
    p = @lift([ Point3f( p_i.p[1],p_i.p[2],p_i.p[3]) for p_i in $cpsO])
    c = @lift([ angle(p_i.p[1]+1im*p_i.p[2]) for p_i in $cpsO])
    arrows3d!(ax, x, p , color=c,  colormap=:phase,colorrange=(-pi,pi), lengthscale=2)
        
    return ax

end

function custom_shape_disks_orientation!(f,ax, cpsO, cfsO)
    shapes = @lift filter(p -> p isa Particles.PolarShapeKernel, $cpsO)
    isempty(shapes[]) && return ax

    c = @lift([ angle(p_i.p[1]+1im*p_i.p[2]) for p_i in $shapes for j in eachindex(p_i.xsrel)])
    s = @lift( [2*p_i.Rs[j] for p_i in $shapes for j in eachindex(p_i.xsrel)] )

    pts = @lift([Point2f(p_i.x[1] + p_i.xsrel[j], p_i.x[2] + p_i.ysrel[j]) for p_i in $shapes for j in eachindex(p_i.xsrel)])

    scatter!(ax, pts, color=c, markersize=s, marker = Circle, markerspace=:data,alpha=1., strokecolor=:transparent, strokewidth=0,colormap=:phase,colorrange=(-pi,pi))
    return ax
end

function simulation()

    R = 5.
    forces = (Forces.self_propulsion(1,0.8),Forces.planar_rotational_noise(ontypes=1,Dr=.0001),
    Forces.shape_repulsive_soft_disk(1,1.),)

    dofevolvers = (DOFevolvers.overdamped_xvf(1),DOFevolvers.overdamped_pqT_xyc(1), DOFevolvers.rigid_xsrel_ysrel(1))

    N=1000
    ϕ = 2

    L =  sqrt(N * pi * R^2 / ϕ)

    r = 1

    xk = [0.]

    #xsrel0 = Float64[ -4/3*sqrt(3)*r ,-1/3*sqrt(3)*r   , 2/3*sqrt(3)*r ,  2/3*sqrt(3)*r , 2/3*sqrt(3)*r , -1/3*sqrt(3)*r]
    #ysrel0 = Float64[ 0, r, 2r, 0, -2r, -r,]

    x_vs = [ -4/3*sqrt(3)*r ,-1/3*sqrt(3)*r   , 2/3*sqrt(3)*r ,  2/3*sqrt(3)*r , 2/3*sqrt(3)*r , -1/3*sqrt(3)*r]
    y_vs = [ 0, r, 2r, 0, -2r, -r,]

    triangle = Initial.construct_populated_polygon(x_vs, y_vs, 2)
    

    square = Initial.construct_populated_polygon([0, 1, 1, 0]*2r, [1, 1, 0, 0]*2r , 2)

    rectangle =  Initial.construct_populated_polygon([0, 2, 2, 0]*2r, [1, 1, 0, 0]*2r , 2)

    shape = square


    xsrel0 = shape.xsrel0 
    rs = r .* ones(length(xsrel0)) 
    display(length(xsrel0))
    #initial_state = ParticleState([ Particles.PolarKernel(id=i,type=1,R =R, x = [rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0],p = normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]), kernel_x = kernel_x, kernel_y = kernel_y) for i=1:N ] );
    initial_state = ParticleState([ Particles.PolarShapeKernel(id=i,type=1, R = R, x = [rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0],p = normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]), xsrel0=shape.xsrel0, ysrel0=shape.ysrel0, Rs=rs) for i=1:N ] );

    system = System(sizes = (L, L, 1.), initial_particle_state = initial_state , forces = forces, dofevolvers=dofevolvers, Periodic = true, rcut_pair_global=2*R);


    sim = Euler_integrator(system,0.01, 1e4,  Tplot=10, res = (1000, 1000),
    plot_functions=( LPlot.shape_orientation!,custom_directors!,), plotdim=2)#

end

sim = simulation()

