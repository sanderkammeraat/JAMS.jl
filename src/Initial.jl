
module Initial
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

end