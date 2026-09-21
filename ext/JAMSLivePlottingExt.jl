module JAMSLivePlottingExt

using JAMS
using GLMakie

function JAMS.GLMakie_window_closeall()
    GLMakie.closeall()
end

function JAMS.setup_system_plotting(system_sizes,plot_functions,plotdim ,cpsO,cfsO,tO,fps;res=nothing,sbs=false)
    GLMakie.activate!(; focus_on_show=true, title= "GLMakie: JAMS simulation", framerate=fps)
    if !isnothing(res)
        f = Figure(size=res)
    else
        f=Figure()
    end
    title = @lift("t = $($tO)")

    if !isnothing(plotdim)
        plotdim_set = plotdim
    else
        plotdim_set = length(system_sizes)

    end

    if plotdim_set==2
        ax = Axis(f[1, 1], xlabel = "x", ylabel="y",  aspect =system_sizes[1]/system_sizes[2], title=title,xgridvisible=false, ygridvisible=false )
        xlims!(ax, -system_sizes[1]/2, system_sizes[1]/2)
        ylims!(ax,  -system_sizes[2]/2, system_sizes[2]/2)

    elseif plotdim_set==3
        if !sbs
        ax = Axis3(f[1, 1], xlabel = "x", ylabel="y", zlabel="z",  aspect = (1,system_sizes[2]/system_sizes[1],system_sizes[3]/system_sizes[1]), title=title)
        xlims!(ax,  -system_sizes[1]/2, system_sizes[1]/2)
        ylims!(ax, -system_sizes[2]/2, system_sizes[2]/2)
        zlims!(ax,  -system_sizes[3]/2, system_sizes[3]/2)

    else
        eye_separation=0.03f0
        ax_left  = Axis3(f[1, 1], xlabel="x", ylabel="y", zlabel="z",aspect=(1,system_sizes[2]/system_sizes[1],system_sizes[3]/system_sizes[1]), title=title)
        ax_right = Axis3(f[1, 2], xlabel="x", ylabel="y", zlabel="z", aspect=(1,system_sizes[2]/system_sizes[1],system_sizes[3]/system_sizes[1]), title=title)
        for ax in (ax_left, ax_right)
                    xlims!(ax, -system_sizes[1]/2, system_sizes[1]/2)
                    ylims!(ax, -system_sizes[2]/2, system_sizes[2]/2)
                    zlims!(ax, -system_sizes[3]/2, system_sizes[3]/2)
        end

        syncing = Ref(false)

        on(ax_left.finallimits) do new_limits
            syncing[] && return
            syncing[] = true
            # Update right axis with the new limits
            xlims!(ax_right, new_limits.origin[1], new_limits.origin[1] + new_limits.widths[1])
            ylims!(ax_right, new_limits.origin[2], new_limits.origin[2] + new_limits.widths[2])
            zlims!(ax_right, new_limits.origin[3], new_limits.origin[3] + new_limits.widths[3])
            syncing[] = false
        end
        
        on(ax_right.finallimits) do new_limits
            syncing[] && return
            syncing[] = true
            # Update left axis with the new limits
            xlims!(ax_left, new_limits.origin[1], new_limits.origin[1] + new_limits.widths[1])
            ylims!(ax_left, new_limits.origin[2], new_limits.origin[2] + new_limits.widths[2])
            zlims!(ax_left, new_limits.origin[3], new_limits.origin[3] + new_limits.widths[3])
            syncing[] = false
        end

        on(ax_left.azimuth) do az
            syncing[] && return
            syncing[] = true
            ax_right.azimuth[] = az + 2f0 * eye_separation  # Reversed sign
            syncing[] = false
        end
        on(ax_right.azimuth) do az
            syncing[] && return
            syncing[] = true
            ax_left.azimuth[] = az - 2f0 * eye_separation  # Reversed sign
            syncing[] = false
        end
        # ── Elevation (identical for both eyes) ───────────────────────────────
        on(ax_left.elevation) do el
            syncing[] && return
            syncing[] = true
            ax_right.elevation[] = el
            syncing[] = false
        end
        on(ax_right.elevation) do el
            syncing[] && return
            syncing[] = true
            ax_left.elevation[] = el
            syncing[] = false
        end
        ax_right.azimuth[] = ax_left.azimuth[] + 2f0 * eye_separation
        end

    end
    if !sbs || plotdim_set==2
        for plot_function in plot_functions
            ax=plot_function(f,ax,cpsO,cfsO)
        end
        display(f)
        return f, ax
    else
        for plot_function in plot_functions
            ax_left=plot_function(f,ax_left,cpsO,cfsO)
            ax_right=plot_function(f,ax_right,cpsO,cfsO)
        end
        display(f)
        return f, (ax_left, ax_right)
    end
end

end

