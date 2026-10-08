
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

#Thanks to Gabriel Martin
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
#Module end
end
