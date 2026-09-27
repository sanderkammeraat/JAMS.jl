
"""
    LPlot

Plot functions for live plotting. Pass a tuple of them to [`Euler_integrator`](@ref JAMS.Euler_integrator) as
`plot_functions`, e.g. `plot_functions=(LPlot.disks_orientation!, LPlot.directors!)`.
Live plotting needs `using GLMakie` after `using JAMS`.

Every plot function has the signature `f(fig, ax, cpsO, cfsO)`, where `cpsO` and `cfsO` are
Makie `Observable`s of the current particle state and field state. Each function only works
for particles (or fields) that have the fields it uses, listed in its docstring.
"""
module LPlot
using CairoMakie
using StaticArrays
using LinearAlgebra
using Distributions
using Random
using SparseArrays

if pkgversion(Makie) >= v"0.23-"
    gen_arrows_2d!(args...; kwargs...) = arrows2d!(args...; kwargs...)
    gen_arrows_3d!(args...; kwargs...) = arrows3d!(args...; kwargs...)
else
    gen_arrows_2d!(args...; kwargs...) = arrows!(args...; kwargs...)
    gen_arrows_3d!(args...; kwargs...) = arrows!(args...; kwargs...)
end


"""
    points!(f, ax, cpsO, cfsO)

Positions, coloured by particle `id`. Requires particle fields `x` and `id`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function points!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ p_i.id[1] for p_i in $cpsO])
    if length(cpsO[][1].x)>2
        
        z = @lift([p_i.x[3] for p_i in $cpsO])
        meshscatter!(ax,x,y,z, color=c)

    else
        scatter!(ax,x,y, color=c)

    end
    return ax
end

#Experimental
"""
    trajectories!(f, ax, cpsO, cfsO)

Experimental. The last 100 positions of every particle, drawn as lines. Requires particle field `x`. Assumes the particle order does not change during the simulation.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function trajectories!(f,ax, cpsO, cfsO)
    maxlen = 100
    N = length(cpsO[])

    trajectories = fill(Point3f(NaN, NaN, NaN), (maxlen+1) * N)
    indices = zeros(Int, N)
    
    obs = Observable(trajectories)
    #Note that particle ids will always be ordered and starting from 1, so we can take 1:N
    color_idx =repeat(1:N, inner=maxlen+1)
    lines!(ax, obs, color=color_idx, colorrange=(1, N), colormap=:viridis)
    
    on(cpsO) do particles
        for i in 1:N
            p = particles[i]
            idx = (indices[i] % maxlen) + 1
            indices[i] = idx
            
            base_idx = (i - 1) * (maxlen+1)
            trajectories[base_idx + idx] = Point3f(p.x[1], p.x[2], p.x[3])
            if idx <maxlen
                trajectories[base_idx + idx+ 1] = Point3f(NaN, NaN, NaN)
            else
                trajectories[base_idx + 2] = Point3f(NaN, NaN, NaN)
            end
    
            
        end
        notify(obs)
    end
    
    return ax
end

"""
    director_points!(f, ax, cpsO, cfsO)

Markers at the tip of the polarity vector, `x + p`, in red (for 3-component positions). Requires particle fields `x`, `p` and `id`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function director_points!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])
    

    c = @lift([ p_i.id[1] for p_i in $cpsO])
    if length(cpsO[][1].x)>2


        
        z = @lift([p_i.x[3] for p_i in $cpsO])
        xpx = @lift([p_i.p[1]+p_i.x[1] for p_i in $cpsO])
        ypy = @lift([p_i.p[2]+p_i.x[2] for p_i in $cpsO])

        meshscatter!(ax,xpx,ypy,z, color=:red)

    else
        scatter!(ax,x,y, color=c)

    end
    return ax
end

"""
    field_magnitude!(f, ax, cpsO, cfsO)

Heatmap of the concentration `C` of the first field, with colour range 0 to 2 and the bin edges drawn as grid lines. Requires a field with fields `bin_centers`, `lbin` and `C`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function field_magnitude!(f,ax, cpsO, cfsO)
    
    field_centers1 = @lift($(cfsO)[1].bin_centers[1])
    field_centers2 = @lift($(cfsO)[1].bin_centers[2])

    field_edges1_l = cfsO[][1].bin_centers[1] .-cfsO[][1].lbin/2
    field_edges2_l = cfsO[][1].bin_centers[2] .-cfsO[][1].lbin/2

    field_edges1_r = cfsO[][1].bin_centers[1] .+cfsO[][1].lbin/2
    field_edges2_r = cfsO[][1].bin_centers[2] .+cfsO[][1].lbin/2
    
    field_C = @lift($(cfsO)[1].C)
    #heatmap!(ax,field_centers1,field_centers2,field_C, alpha=0.2,colormap=:jet,colorrange=(0.0,1))
    #Colorbar(f[1,2], limits = (0.0, 1), label="Concentration c",colormap=(:jet,0.2))

    heatmap!(ax,field_centers1,field_centers2,field_C, alpha=0.2,colormap=:jet,colorrange=(0.0,2))
    Colorbar(f[1,2], limits = (0, 2), label="Concentration c",colormap=(:jet,0.2))
    vlines!(ax,field_edges1_l, color="grey", alpha=0.1)
    hlines!(ax,field_edges2_l, color="grey", alpha=0.1)
    vlines!(ax,field_edges1_r, color="grey", alpha=0.1)
    hlines!(ax,field_edges2_r, color="grey", alpha=0.1)
    return ax
end

"""
    GPUfield_magnitude!(f, ax, cpsO, cfsO)

Same as [`field_magnitude!`](@ref), but converts `C` to `Float64` first, for fields stored in single precision or on a GPU.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function GPUfield_magnitude!(f,ax, cpsO, cfsO)
    
    field_centers1 = @lift($(cfsO)[1].bin_centers[1])
    field_centers2 = @lift($(cfsO)[1].bin_centers[2])

    field_edges1_l = cfsO[][1].bin_centers[1] .-cfsO[][1].lbin/2
    field_edges2_l = cfsO[][1].bin_centers[2] .-cfsO[][1].lbin/2

    field_edges1_r = cfsO[][1].bin_centers[1] .+cfsO[][1].lbin/2
    field_edges2_r = cfsO[][1].bin_centers[2] .+cfsO[][1].lbin/2
    
    field_C = @lift(Float64.($(cfsO)[1].C))
    #heatmap!(ax,field_centers1,field_centers2,field_C, alpha=0.2,colormap=:jet,colorrange=(0.0,1))
    #Colorbar(f[1,2], limits = (0.0, 1), label="Concentration c",colormap=(:jet,0.2))

    heatmap!(ax,field_centers1,field_centers2,field_C, alpha=0.2,colormap=:jet,colorrange=(0.0,2))
    Colorbar(f[1,2], limits = (0, 2), label="Concentration c",colormap=(:jet,0.2))
    vlines!(ax,field_edges1_l, color="grey", alpha=0.1)
    hlines!(ax,field_edges2_l, color="grey", alpha=0.1)
    vlines!(ax,field_edges1_r, color="grey", alpha=0.1)
    hlines!(ax,field_edges2_r, color="grey", alpha=0.1)
    return ax
end

"""
    field_log_magnitude!(f, ax, cpsO, cfsO)

Heatmap of `log10(C)` of the first field, with colour range -4 to -2. Requires a field with fields `bin_centers` and `C`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function field_log_magnitude!(f,ax, cpsO, cfsO)
    
    field_centers1 = @lift($(cfsO)[1].bin_centers[1])
    field_centers2 = @lift($(cfsO)[1].bin_centers[2])


    
    field_C = @lift(log10.($(cfsO)[1].C))
    heatmap!(ax,field_centers1,field_centers2,field_C, alpha=0.2,colormap=:jet,colorrange=(-4,-2))
    Colorbar(f[1,2], limits = (-4, -2), label="Concentration c",colormap=(:jet,0.2))
    return ax
end



"""
    potential!(f, ax, cpsO, cfsO)

Surface plot of `C` of the first field, interpreted as a potential energy landscape. Use with `plotdim=3`. Requires a field with fields `bin_centers` and `C`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function potential!(f,ax, cpsO, cfsO)
    
    field_centers1 = @lift($(cfsO)[1].bin_centers[1])
    field_centers2 = @lift($(cfsO)[1].bin_centers[2])
    field_C = @lift($(cfsO)[1].C)
    Cmax = maximum(cfsO[][1].C)
    Cmin = minimum(cfsO[][1].C)
    surface!(ax,field_centers1,field_centers2,field_C,colormap=Reverse(:gist_rainbow), alpha=0.5, colorrange=(Cmin, Cmax), color=field_C,shading=NoShading)
    Colorbar(f[1,2], limits = (Cmin, Cmax), label="Pot. energy U",colormap=Reverse(:gist_rainbow))
    return ax
end

"""
    field_magnitude_wgrid!(f, ax, cpsO, cfsO)

Heatmap of `C` of the first field, with colour range 0 to 1 and grid lines at the bin centres. Requires a field with fields `bin_centers` and `C`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function field_magnitude_wgrid!(f,ax, cpsO, cfsO)
    
    field_centers1 = @lift($(cfsO)[1].bin_centers[1])
    field_centers2 = @lift($(cfsO)[1].bin_centers[2])
    field_C = @lift($(cfsO)[1].C)
    heatmap!(ax,field_centers1,field_centers2,field_C, alpha=0.2,colormap=:viridis,colorrange=(0,1))

    vlines!(ax,field_centers1, color="black", alpha=0.2)
    hlines!(ax,field_centers2, color="black", alpha=0.2)
    return ax
end

"""
    type_points!(f, ax, cpsO, cfsO)

Positions, coloured by particle `type`. Requires particle fields `x` and `type`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function type_points!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ p_i.type[1] for p_i in $cpsO])
    if length(cpsO[][1].x)>2
        
        z = @lift([p_i.x[3] for p_i in $cpsO])
        meshscatter!(ax,x,y,z, color=c)

    else
        scatter!(ax,x,y, color=c)

    end
    return ax
end


"""
    shape_disks!(f, ax, cpsO, cfsO)

Rigid bodies drawn as circles at their extent points `xe` with radii `re`, coloured by `id`. Requires particle fields `xe`, `re` and `id`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function shape_disks!(f,ax, cpsO, cfsO)
    c = @lift([ p_i.id[1] for p_i in $cpsO])
    for j = 1:size(cpsO[][1].xe)[1]

        xes = @lift([p_i.xe[j,1] for p_i in $cpsO] )

        yes = @lift([p_i.xe[j,2] for p_i in $cpsO] )

        s = @lift([2*p_i.re[j] for p_i in $cpsO])
    
        scatter!(ax,xes,yes, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.5, strokecolor=:black, strokewidth=1)
    end
    return ax
end



"""
    shape_disks_orientation!(f, ax, cpsO, cfsO)

Like [`shape_disks!`](@ref), coloured by the angle of the polarity `p` in the xy-plane.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function shape_disks_orientation!(f,ax, cpsO, cfsO)
    c = @lift([ angle(p_i.p[1]+1im*p_i.p[2]) for p_i in $cpsO])
    for j = 1:size(cpsO[][1].xe)[1]

        xes = @lift([p_i.xe[j,1] for p_i in $cpsO] )

        yes = @lift([p_i.xe[j,2] for p_i in $cpsO] )

        s = @lift([2*p_i.re[j] for p_i in $cpsO])
    
        scatter!(ax,xes,yes, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.5, strokecolor=:black, strokewidth=1,colormap=:hsv,colorrange=(-pi,pi))
    end
    return ax
end

"""
    shape_disks_type!(f, ax, cpsO, cfsO)

Like [`shape_disks!`](@ref), coloured by particle `type`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function shape_disks_type!(f,ax, cpsO, cfsO)
    c = @lift([ p_i.type[1] for p_i in $cpsO])
    for j = 1:size(cpsO[][1].xe)[1]

        xes = @lift([p_i.xe[j,1] for p_i in $cpsO] )

        yes = @lift([p_i.xe[j,2] for p_i in $cpsO] )

        s = @lift([2*p_i.re[j] for p_i in $cpsO])
    
        scatter!(ax,xes,yes, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.5, strokecolor=:black, strokewidth=1,colormap=Reverse(:seismic))
    end
    return ax
end

"""
    shape_points!(f, ax, cpsO, cfsO)

Rigid bodies drawn as points at their extent points `xe`, coloured by `id`. Requires particle fields `xe`, `re` and `id`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function shape_points!(f,ax, cpsO, cfsO)
    c = @lift([ p_i.id[1] for p_i in $cpsO])
    for j = 1:size(cpsO[][1].xe)[1]

        xes = @lift([p_i.xe[j,1] for p_i in $cpsO] )

        yes = @lift([p_i.xe[j,2] for p_i in $cpsO] )

        s = @lift([p_i.re[1] for p_i in $cpsO])
    
        scatter!(ax,xes,yes, color=c)
    end
    return ax
end




"""
    Swarmalators!(f, ax, cpsO, cfsO)

Positions, coloured by the internal phase `ϕ` (swarmalator models). Requires particle fields `x` and `ϕ`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function Swarmalators!(f,ax, cpsO, cfsO)

    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift(angle2range.([ p_i.ϕ[1] for p_i in $cpsO]))

    scatter!(ax,x,y , color=c, colormap=:hsv, colorrange=(0, 2*pi))
    return ax
end


"""
    sphere!(f, ax, cpsO, cfsO)

A white sphere centred on the origin, with radius equal to the distance of the first particle from the origin. Useful as a background for particles confined to a sphere. Use with `plotdim=3`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function sphere!(f,ax, cpsO, cfsO)


    R = norm(cpsO[][1].x)
    meshscatter!(ax,0,0,0, markersize=R, color=:white)

    return ax
end


"""
    sized_points!(f, ax, cpsO, cfsO)

Particles drawn at their true size: disks of radius `R` in 2D, transparent spheres in 3D, coloured by `id`. Requires particle fields `x`, `R` and `id`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function sized_points!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ p_i.id[1] for p_i in $cpsO])

    

    if length(cpsO[][1].x)>2
        
        z = @lift([p_i.x[3] for p_i in $cpsO])

        R = @lift([p_i.R[1] for p_i in $cpsO])

        meshscatter!(ax,x,y,z, color=c, markersize =R, transparency=true, colormap=(:viridis,0.2))


    else
        s = @lift([2*p_i.R[1]  for p_i in $cpsO])
        scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1)

    end

    return ax
end


"""
    polymers!(f, ax, cpsO, cfsO)

Disks of radius `R`, coloured by polymer id `pol_id`. Requires particle fields `x`, `R` and `pol_id`, as in [`Particles.PolarPolymer`](@ref JAMS.Particles.PolarPolymer).

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function polymers!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ p_i.pol_id[1] for p_i in $cpsO])

    
    s = @lift([2*p_i.R[1]  for p_i in $cpsO])
    scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1, colormap=:rainbow1)

    return ax
end


"""
    polymers_3d!(f, ax, cpsO, cfsO)

Spheres of radius `R`, coloured by polymer id `pol_id`. Use with `plotdim=3`. Requires particle fields `x`, `R` and `pol_id`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function polymers_3d!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])
    z = @lift([p_i.x[3] for p_i in $cpsO])

    c = @lift([ p_i.pol_id[1] for p_i in $cpsO])

    
    s = @lift([p_i.R[1]  for p_i in $cpsO])
    meshscatter!(ax,x,y,z, color=c, markersize =s,alpha=0.7, colormap=:rainbow1)

    return ax
end

"""
    disks!(f, ax, cpsO, cfsO)

Disks of radius `R` at the particle positions, coloured by `id`. Requires particle fields `x`, `R` and `id`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function disks!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ p_i.id[1] for p_i in $cpsO])

    
    s = @lift([2*p_i.R[1]  for p_i in $cpsO])
    scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1)

    return ax
end



"""
    ellipse(p_i) -> (xs, ys)

Outline of the ellipse used by [`ellipses!`](@ref) for particle `p_i`.
"""
function ellipse(p_i)

    lmda_major = 0.5*(p_i.Lambda[1,1]+p_i.Lambda[2,2]) + sqrt(0.25*(p_i.Lambda[1,1]-p_i.Lambda[2,2])^2 + p_i.Lambda[1,2]^2)

    lmda_minor =0.5*(p_i.Lambda[1,1]+p_i.Lambda[2,2]) -sqrt(0.25*(p_i.Lambda[1,1]-p_i.Lambda[2,2])^2 + p_i.Lambda[1,2]^2)

    s = range(0, 2pi, length=40)

    x_e = lmda_major .* cos.(s)

    y_e = lmda_minor .* sin.(s)

   x_e_rot = p_i.p[1] * x_e  - p_i.p[2] * y_e
   y_e_rot= p_i.p[2] * x_e  + p_i.p[1] *  y_e

   return p_i.x[1] .+ x_e_rot, p_i.x[2] .+ y_e_rot

end



"""
    ellipses!(f, ax, cpsO, cfsO)

Particles drawn as ellipses whose axes are the eigenvalues of the 2×2 shape tensor `Lambda`, rotated along the polarity `p`. Requires particle fields `x`, `p` and `Lambda`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function ellipses!(f,ax, cpsO, cfsO)

    xs = @lift([ ellipse( p_i )[1] for p_i in $cpsO])
    ys = @lift([ ellipse( p_i )[2] for p_i in $cpsO])


    for i=1:length(cpsO[])

        xsi = @lift( $xs[i])
        ysi = @lift($ys[i])

        poly!(ax, xsi, ysi, color=i, colorrange=(1,length(cpsO[])),alpha=0.7, strokecolor=:black, strokewidth=2)
    end

    return ax
end

"""
    transparant_disks!(f, ax, cpsO, cfsO)

White, almost transparent disks of radius `R` with a thin outline. Requires particle fields `x` and `R`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function transparant_disks!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])
    
    s = @lift([2*p_i.R[1]  for p_i in $cpsO])
    scatter!(ax,x,y, color=:white, markersize =s,marker = Circle, markerspace=:data,alpha=0.1, strokecolor=:black, strokewidth=.2)

    return ax
end



"""
    disks_type!(f, ax, cpsO, cfsO)

Disks of radius `R`, coloured by particle `type`. Requires particle fields `x`, `R` and `type`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function disks_type!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ p_i.type[1] for p_i in $cpsO])

    
    s = @lift([2*p_i.R[1]  for p_i in $cpsO])
    scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1, colormap=:RdYlGn)

    return ax
end


"""
    disks_uw!(f, ax, cpsO, cfsO)

Disks of radius `R` at the *unwrapped* positions `xuw`, coloured by `id`. Requires particle fields `xuw`, `R` and `id`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function disks_uw!(f,ax, cpsO, cfsO)


    x = @lift([p_i.xuw[1] for p_i in $cpsO])
    y = @lift([p_i.xuw[2] for p_i in $cpsO])

    c = @lift([ p_i.id[1] for p_i in $cpsO])

    
    s = @lift([2*p_i.R[1]  for p_i in $cpsO])
    scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1)

    return ax
end

"""
    disks_vx!(f, ax, cpsO, cfsO)

Disks of radius `R`, coloured by the x-component of the velocity (colour range ±0.01). Requires particle fields `x`, `R` and `v`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function disks_vx!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ p_i.v[1] for p_i in $cpsO])

    
    s = @lift([2*p_i.R[1]  for p_i in $cpsO])
    scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1,colorrange=(-0.01,0.01),colormap=:seismic)

    return ax
end


"""
    disks_orientation!(f, ax, cpsO, cfsO)

Disks of radius `R`, coloured by the angle of the polarity `p` in the xy-plane. Requires particle fields `x`, `R` and `p`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function disks_orientation!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ angle(p_i.p[1]+1im*p_i.p[2]) for p_i in $cpsO])

    
    s = @lift([2*p_i.R[1]  for p_i in $cpsO])
    scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1,colormap=:hsv,colorrange=(-pi,pi))

    return ax
end

"""
    disks_v_orientation!(f, ax, cpsO, cfsO)

Disks of radius `R`, coloured by the direction of the velocity `v` in the xy-plane. Requires particle fields `x`, `R` and `v`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function disks_v_orientation!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ angle(p_i.v[1]+1im*p_i.v[2]) for p_i in $cpsO])

    
    s = @lift([2*p_i.R[1]  for p_i in $cpsO])
    scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1,colormap=:hsv,colorrange=(-pi,pi))

    return ax
end

"""
    disks_nematic_orientation!(f, ax, cpsO, cfsO)

Disks of radius `R`, coloured by the nematic angle of the polarity (the angle of `p` doubled, so `p` and `-p` get the same colour). Requires particle fields `x`, `R` and `p`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function disks_nematic_orientation!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ angle( exp(angle(p_i.p[1]+1im*p_i.p[2])*2  * 1im)) for p_i in $cpsO])

    
    s = @lift([2*p_i.R[1]  for p_i in $cpsO])
    scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1,colormap=:hsv,colorrange=(-pi,pi))

    return ax
end

"""
    disks_vp_phase_difference!(f, ax, cpsO, cfsO)

Disks of radius `R`, coloured by the angle between the velocity `v` and the polarity `p`. Requires particle fields `x`, `R`, `v` and `p`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function disks_vp_phase_difference!(f,ax, cpsO, cfsO)


    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ angle( exp(1im* (angle(p_i.v[1]+1im*p_i.v[2]) - angle(p_i.p[1]+1im*p_i.p[2])) ) )  for p_i in $cpsO])

    
    s = @lift([2*p_i.R[1]  for p_i in $cpsO])
    scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1,colormap=:hsv,colorrange=(-pi,pi))

    return ax
end


"""
    type_sized_points!(f, ax, cpsO, cfsO)

Particles drawn at their true size (disks in 2D, spheres in 3D), coloured by `type`. Requires particle fields `x`, `R` and `type`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function type_sized_points!(f,ax, cpsO, cfsO)

    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    c = @lift([ p_i.type[1] for p_i in $cpsO])

    

    if length(cpsO[][1].x)>2
        
        z = @lift([p_i.x[3] for p_i in $cpsO])

        R = @lift([p_i.R for p_i in $cpsO])

        meshscatter!(ax,x,y,z, color=c, markersize =R, transparency=true,markerspace=:data)


    else
        s = @lift([2*p_i.R[1]  for p_i in $cpsO])
        scatter!(ax,x,y, color=c, markersize =s,marker = Circle, markerspace=:data,alpha=0.7, strokecolor=:black, strokewidth=1)

    end
    return ax

end

"""
    directors!(f, ax, cpsO, cfsO)

Arrows along the polarity `p` at every particle, coloured by the angle of `p` in the xy-plane. Requires particle fields `x` and `p`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function directors!(f,ax, cpsO, cfsO)

    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    
    if length(cpsO[][1].x)>2

        x = @lift([ Point3f( p_i.x[1],p_i.x[2],p_i.x[3]) for p_i in $cpsO])
        p = @lift([ Point3f( p_i.p[1],p_i.p[2],p_i.p[3]) for p_i in $cpsO])
        c = @lift([ angle(p_i.p[1]+1im*p_i.p[2]) for p_i in $cpsO])
        gen_arrows_3d!(ax, x, p , color=c,  colormap=:hsv,colorrange=(-pi,pi))
        


    else
        nx = @lift(cos.([ p_i.θ[1] for p_i in $cpsO]))
        ny = @lift(sin.([ p_i.θ[1] for p_i in $cpsO]))
        c = @lift(angle2range.([ p_i.θ[1] for p_i in $cpsO]))
        gen_arrows_2d!(ax, x, y, nx, ny, color=c,  colormap=:hsv,colorrange=(0, 2*pi))


    end
    return ax

end

"""
    nematic_directors!(f, ax, cpsO, cfsO)

Double-headed arrows along `p` and `-p` at every particle, coloured by the nematic angle. Requires particle fields `x` and `p`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function nematic_directors!(f,ax, cpsO, cfsO)

    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    
    if length(cpsO[][1].x)>2

        x = @lift([ Point3f( p_i.x[1],p_i.x[2],p_i.x[3]) for p_i in $cpsO])
        p = @lift([ Point3f( p_i.p[1],p_i.p[2],p_i.p[3]) for p_i in $cpsO])

        pmin = @lift(-1*[ Point3f( p_i.p[1],p_i.p[2],p_i.p[3]) for p_i in $cpsO])


        c = @lift([ angle( exp(angle(p_i.p[1]+1im*p_i.p[2])*2  * 1im)) for p_i in $cpsO])
        gen_arrows_3d!(ax, x, p , color=c,  colormap=:hsv,colorrange=(-pi,pi))
        gen_arrows_3d!(ax, x, pmin , color=c,  colormap=:hsv,colorrange=(-pi,pi))


    else
        nx = @lift(cos.([ p_i.θ[1] for p_i in $cpsO]))
        ny = @lift(sin.([ p_i.θ[1] for p_i in $cpsO]))


        nxmin = @lift(-1*cos.([ p_i.θ[1] for p_i in $cpsO]))
        nymin = @lift(-1*sin.([ p_i.θ[1] for p_i in $cpsO]))


        c = @lift([ angle( exp(angle(p_i.p[1]+1im*p_i.p[2])*2  * 1im)) for p_i in $cpsO])
        gen_arrows_2d!(ax, x, y, nx, ny, color=c,  colormap=:hsv,colorrange=(0, 2*pi))
        gen_arrows_2d!(ax, x, y, nxmin, nymin, color=c,  colormap=:hsv,colorrange=(0, 2*pi))

    end
    return ax

end

"""
    velocity_vectors!(f, ax, cpsO, cfsO)

Arrows along the velocity `v` at every particle, coloured by the direction of `v`. Requires particle fields `x` and `v`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function velocity_vectors!(f,ax,cpsO, cfsO)

    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    
    if length(cpsO[][1].x)>2

        x = @lift([ Point3f( p_i.x[1],p_i.x[2],p_i.x[3]) for p_i in $cpsO])
        v = @lift([ Point3f( p_i.v[1],p_i.v[2],p_i.v[3]) for p_i in $cpsO])
        c = @lift([ angle(p_i.v[1]+1im*p_i.v[2]) for p_i in $cpsO])
        gen_arrows_3d!(ax, x, v , color=c,  colormap=:hsv,colorrange=(-pi,pi))
        


    else
        vx = @lift([ p_i.v[1] for p_i in $cpsO])
        vy = @lift([ p_i.v[2] for p_i in $cpsO])
        c = @lift([ angle(p_i.v[1]+1im*p_i.v[2]) for p_i in $cpsO])
        gen_arrows_2d!(ax, x, y, vx, vy, color=c, colormap=:hsv,colorrange=(-pi,pi))


    end
    return ax
end

"""
    annotate_v!(f, ax, cpsO, cfsO)

The speed `|v|` of every particle as a text label next to it. Requires particle fields `x` and `v`.

Live plot function: pass it to [`Euler_integrator`](@ref JAMS.Euler_integrator) in `plot_functions`.
"""
function annotate_v!(f,ax,cpsO, cfsO)

    x = @lift([p_i.x[1] for p_i in $cpsO])
    y = @lift([p_i.x[2] for p_i in $cpsO])

    
    vtext = @lift([ string(round( sqrt( p_i.v[1]^2 + p_i.v[2]^2 + p_i.v[3]^2); digits=2 )) for p_i in $cpsO])

    text!(ax, x,y, text=vtext)




    return ax
end




"""
    angle2range(angle)

Map an angle in radians to the range `[0, 2π)`.
"""
function angle2range(angle)
    if angle>=0
        return angle % (2*pi)
    else
        return angle % (2*pi) + (2*pi)
    end

end



# function setup_system_plotting(system_sizes,plot_functions,plotdim ,cpsO,cfsO,tO,fps;res=nothing)
#     GLMakie.activate!(; focus_on_show=true, title= "GLMakie: JAMs simulation", framerate=fps)
#     if !isnothing(res)
#         f = Figure(size=res)
#     else
#         f=Figure()
#     end
#     title = @lift("t = $($tO)")

#     if !isnothing(plotdim)
#         plotdim_set = plotdim
#     else
#         plotdim_set = length(system_sizes)

#     end

#     if plotdim_set==2
#         ax = Axis(f[1, 1], xlabel = "x", ylabel="y",  aspect =system_sizes[1]/system_sizes[2], title=title )
#         xlims!(ax, -system_sizes[1]/2, system_sizes[1]/2)
#         ylims!(ax,  -system_sizes[2]/2, system_sizes[2]/2)

#     elseif plotdim_set==3

#         ax = Axis3(f[1, 1], xlabel = "x", ylabel="y", zlabel="z",  aspect = (1,system_sizes[2]/system_sizes[1],system_sizes[3]/system_sizes[1]), title=title)
#         xlims!(ax,  -system_sizes[1]/2, system_sizes[1]/2)
#         ylims!(ax, -system_sizes[2]/2, system_sizes[2]/2)
#         zlims!(ax,  -system_sizes[3]/2, system_sizes[3]/2)

#     end
    
#     for plot_function in plot_functions
#        ax=plot_function(f,ax,cpsO,cfsO)
#     end
#     display(f)
#     return f, ax
# end


end