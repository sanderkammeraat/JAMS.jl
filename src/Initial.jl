
module Initial
using StaticArrays
using LinearAlgebra
using Distributions
using Random
using SparseArrays
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


function random_in_disk(N, R_out;R_in=0)


    rs = rand(Uniform(0,1),N)

    thetas = rand(Uniform(0,2pi),N)

    xs = sqrt.(rs .* (R_out^2 - R_in^2 ) .+ R_in^2 ) .* cos.(thetas)
    ys = sqrt.(rs .* (R_out^2 - R_in^2 ) .+ R_in^2 ) .* sin.(thetas)

    return xs, ys

end

#Module end
end