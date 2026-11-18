
module Initial
using StaticArrays
using LinearAlgebra
using Distributions
using Random
using SparseArrays

"""
    box(l, Lx, Ly) -> (xs, ys)

Positions spaced roughly `l` apart along the edge of a rectangle of size `Lx × Ly`
centred on the origin.

Returns two vectors with the x- and y-coordinates.
"""
function box(l, Lx, Ly)

    x_range = collect(range(-Lx/2, Lx/2,length=trunc(Int64,cld(Lx,l))))
    Nx = length(x_range)

    y_range = collect(range(-Ly/2, Ly/2,length=trunc(Int64,cld(Ly,l))))
    Ny = length(y_range)

    xs = []
    ys = []

    for i=1:Nx-1
        push!(xs, x_range[i])
        push!(ys, y_range[end])
    end
    for i=1:Ny-1
        push!(xs, x_range[end])
        push!(ys, y_range[end-i+1])
    end
    for i=1:Nx-1
        push!(xs, x_range[end-i+1])
        push!(ys, y_range[1])
    end
    for i=1:Ny-1
        push!(xs, x_range[1])
        push!(ys, y_range[i])
    end

    return xs, ys

end

"""
    ring(l, R) -> (xs, ys)

Positions on a circle of radius `R` centred on the origin, spaced an arc length of about
`l` apart. Returns empty vectors if `R == 0`.

Returns two vectors with the x- and y-coordinates.
"""
function ring(l,R)

    if R!=0
    dtheta = l/R

    thetas = 0:dtheta:2pi

    xs = R .*cos.(thetas)
    ys = R .*sin.(thetas)
    else
        xs = []
        ys = []
    end
    return xs, ys

end


"""
    random_in_disk(N, R_out; R_in=0) -> (xs, ys)

`N` random positions distributed uniformly (by area) over a disk of radius `R_out`
centred on the origin, or over the annulus `R_in ≤ r ≤ R_out` if `R_in > 0`.

Returns two vectors with the x- and y-coordinates.
"""
function random_in_disk(N, R_out;R_in=0)


    rs = rand(Uniform(0,1),N)

    thetas = rand(Uniform(0,2pi),N)

    xs = sqrt.(rs .* (R_out^2 - R_in^2 ) .+ R_in^2 ) .* cos.(thetas)
    ys = sqrt.(rs .* (R_out^2 - R_in^2 ) .+ R_in^2 ) .* sin.(thetas)

    return xs, ys

end


"""
    stacked_polymers_at_angle(N_in_pol, R, pf, f_eq_stretch_force, L_0; tilt_angle=nothing, random_polarity=false)
        -> (x, y, radii, pol_ids, ids_in_pol, L, Npols)

Credit: G. Martin. 
Straight polymers of `N_in_pol` monomers of radius `R` in a periodic square box, at
packing fraction `pf`. Consecutive monomers are `2R * f_eq_stretch_force` apart (use the
same value as `farray` in [`Forces.polymer_harmonic_stretch`](@ref JAMS.Forces.polymer_harmonic_stretch)),
and the overlap between them is taken into account in the packing fraction.

The number of polymers `Npols` is chosen so that the box side `L` is close to `L_0`. The
polymers are placed end to end along one long line at angle `tilt_angle` to the x-axis,
wrapped around the periodic box. By default `tilt_angle` is chosen so that the wrapped line
spreads the polymers evenly over the box. Polymer ids are assigned in random order. If
`random_polarity` is `true`, the monomer order of each polymer is reversed with probability
1/2, so the polarity of those polymers (see [`DOFevolvers.polymer_p_set`](@ref JAMS.DOFevolvers.polymer_p_set))
points the other way.

Returns the x- and y-coordinates, radii, polymer ids and positions along the polymer of all
`Npols * N_in_pol` monomers, plus the box side `L` and `Npols`. Use `L` for the box size
and set `Periodic = true`:

```julia
x, y, radii, pol_ids, ids_in_pol, L, Npols = Initial.stacked_polymers_at_angle(10, 1.0, 0.9, 0.75, 80.0)
initial_state = ParticleState([Particles.PolarPolymer(id=i, type=1, pol_id=pol_ids[i],
    id_in_pol=ids_in_pol[i], pol_N=10, R=radii[i], x=[x[i], y[i], 0.0], p=[1.0, 0.0, 0.0])
    for i in 1:Npols*10])
```
"""
function stacked_polymers_at_angle(N_in_pol, R, pf, f_eq_stretch_force, L_0; tilt_angle = nothing, random_polarity = false)

    # initialization

    S = pi * R^2

    S_overlap = 2*(R*R*acos(f_eq_stretch_force*R) - (f_eq_stretch_force*R)^2*tan(acos(f_eq_stretch_force*R)))

    S_polymers = (S*N_in_pol-S_overlap*(N_in_pol-1))

    Npols = floor(Int64, pf * L_0^2 / S_polymers)

    L = sqrt(Npols*S_polymers/pf)

    M = Npols*N_in_pol

    x = zeros(M)
    y = zeros(M)
    radii = R * ones(M)

    pol_ids = zeros(M)

    ids_in_pol= zeros(M)

    Lpols = 2*R + (N_in_pol - 1) * 2*R * f_eq_stretch_force

    if isnothing(tilt_angle)
        tilt_angle = atan(L/Npols/Lpols)   # The angle along which to place the particles to maximize the spread of particles on the torque
    end
    
    id = 1

    for (i, pol_index) in enumerate(shuffle!(collect(1:Npols)))

        flip = random_polarity ? rand((true,false)) : false   #based on the random_polarity bool, if true choose random, if false use false

        for j in 1:N_in_pol
            
            pol_ids[id] = pol_index

            ids_in_pol[id] = j

            index = j - 1

            if flip
                index = N_in_pol - j
            end

            x[id] = ((i-1)*Lpols + index*2*R*f_eq_stretch_force)*cos(tilt_angle) % L - L/2
            y[id] = ((i-1)*Lpols + index*2*R*f_eq_stretch_force)*sin(tilt_angle) % L - L/2

            id+=1

        end

    end
    return x, y, radii, pol_ids, ids_in_pol, L, Npols

end


function circular_cap(opening_angle, R,lbin)

    Nx = ceil(Int64, 2*R/lbin)
    Ny = copy(Nx)

    x = collect(0:Nx-1) .* lbin .- R .+ lbin/2
    y = collect(0:Ny-1) .* lbin .- R .+ lbin/2

    X = zeros(Nx, Ny)

    Y = zeros(Nx, Ny)

    square = zeros(Bool, Nx, Ny)

    for (i, xi)  in pairs(x)
        for (j,yj) in pairs(y)

            X[i,j] = x[i]

            Y[i,j] = y[j]
        end
    end

    angles = atan.(Y,X)
    
    square[ (X.^2 + Y.^2 .< (R + lbin/2)^2 )   .*  (X.^2 + Y.^2 .> (R - 1*lbin).^2)   .* ( (angles .< -pi+opening_angle) .| (angles .> pi-opening_angle)  )] .= true

    xk = X[square]
    yk = Y[square]

    return xk, yk;

end


#Polygon constructions using edges
@kwdef struct Edge
    x1::Float64
    y1::Float64

    x2::Float64
    y2::Float64

    Lx = x2 - x1
    Ly = y2 - y1

    L = sqrt(Lx^2 + Ly^2)
end

struct PopulatedEdge
    xs
    ys
end

function place_on_edge(Edge, density)

    N = ceil(Int, Edge.L * density)

    dLx = Edge.Lx/(N-1)
    dLy = Edge.Ly/(N-1)

    xs = [Edge.x1 + dLx*n for n=0:N]

    ys = [Edge.y1 + dLy*n for n=0:N]

    return PopulatedEdge(xs, ys)

end

struct PopulatedPolygon

    xcom
    ycom

    xs
    ys

    xsrel0
    ysrel0

end

function construct_populated_polygon_from_edges(populated_edges)

    xp = []
    yp = []

    for (i, populated_edge) in enumerate(populated_edges)

        #Loop clockwise, discard the last point of each populated edge, to avoid duplicates
        for j in 1:length(populated_edge.xs) - 1

            push!(xp, populated_edge.xs[j])
            push!(yp, populated_edge.ys[j])
        end


    end

    xcom = mean(xp)
    ycom = mean(yp)

    xprel = xp .- xcom
    yprel = yp .- ycom

    return PopulatedPolygon(xcom, ycom, xp, yp, xprel, yprel)

end

#Unique coordinates of the vertices only, ordered counterclockwise
function construct_populated_polygon(x_vertices, y_vertices, density)

    Nv = length(x_vertices)

    edges = [ Edge(x1 = x_vertices[i], y1 = y_vertices[i], x2 = x_vertices[i+1], y2 = y_vertices[i+1] ) for i=1:Nv-1]

    push!(edges, Edge(x1 = x_vertices[Nv], y1 = y_vertices[Nv], x2 = x_vertices[1], y2 = y_vertices[1]))

    populated_edges = [ place_on_edge(edges[i], density) for i in eachindex(edges)]


    populated_polygon = construct_populated_polygon_from_edges(populated_edges)

    return populated_polygon
end







#Module end
end