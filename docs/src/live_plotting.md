# Live plotting

JAMS can show the simulation while it runs, using GLMakie. This is useful for quickly
seeing how a model behaves for different parameters. Load GLMakie *after* JAMS to enable
it:

```julia
using JAMS
using GLMakie
```

Then give [`Euler_integrator`](@ref) a plot interval and a tuple of plot functions:

```julia
sim = Euler_integrator(system, 0.01, 1000.0;
    Tplot = 10,                                                  # redraw every 10 timesteps
    fps = 60,
    plot_functions = (LPlot.disks_orientation!, LPlot.directors!),
)
```

The plot functions are drawn on top of each other in the order given, so here the
particles are drawn as disks coloured by orientation, with arrows along their polarity on
top. Closing the window stops the simulation.

Other plotting keywords of [`Euler_integrator`](@ref):

- `plotdim = 3`: use a 3D axis, for 3D simulations (**Default**: `2`)
- `sbs = true`: with `plotdim = 3`, show two side-by-side 3D views for stereo viewing
- `res = (1000, 1000)`: window size in pixels
- `record_folder_path`: also record the plot to `movie.mp4` in this folder; set the format
  with `format` and the compression with `crf`

On a cluster without a display, leave out `using GLMakie` and `Tplot`.

Each plot function only works for particles (or fields) that have the fields it uses. The
tables below list them, and the [`LPlot`](@ref) module docstring describes the signature
for writing your own; see also [A custom plot function](@ref).

## Particles as points

| Plot function | Shows | Needs particle fields |
|---|---|---|
| [`LPlot.points!`](@ref) | positions, coloured by id | `x`, `id` |
| [`LPlot.type_points!`](@ref) | positions, coloured by type | `x`, `type` |
| [`LPlot.director_points!`](@ref) | markers at the tip of the polarity, `x + p` | `x`, `p`, `id` |
| [`LPlot.trajectories!`](@ref) | the last 100 positions of every particle as lines (experimental) | `x` |
| [`LPlot.Swarmalators!`](@ref) | positions, coloured by internal phase | `x`, `ϕ` |

## Particles at their true size

| Plot function | Shows | Needs particle fields |
|---|---|---|
| [`LPlot.disks!`](@ref) | disks of radius `R`, coloured by id | `x`, `R`, `id` |
| [`LPlot.disks_type!`](@ref) | disks, coloured by type | `x`, `R`, `type` |
| [`LPlot.disks_orientation!`](@ref) | disks, coloured by polarity angle | `x`, `R`, `p` |
| [`LPlot.disks_nematic_orientation!`](@ref) | disks, coloured by nematic angle (`p` and `-p` same colour) | `x`, `R`, `p` |
| [`LPlot.disks_v_orientation!`](@ref) | disks, coloured by velocity direction | `x`, `R`, `v` |
| [`LPlot.disks_vp_phase_difference!`](@ref) | disks, coloured by the angle between `v` and `p` | `x`, `R`, `v`, `p` |
| [`LPlot.disks_vx!`](@ref) | disks, coloured by the x-component of the velocity | `x`, `R`, `v` |
| [`LPlot.disks_uw!`](@ref) | disks at the unwrapped positions | `xuw`, `R`, `id` |
| [`LPlot.transparant_disks!`](@ref) | white, almost transparent disks with an outline | `x`, `R` |
| [`LPlot.sized_points!`](@ref) | disks in 2D, transparent spheres in 3D, coloured by id | `x`, `R`, `id` |
| [`LPlot.type_sized_points!`](@ref) | disks in 2D, spheres in 3D, coloured by type | `x`, `R`, `type` |
| [`LPlot.polymers!`](@ref) | disks, coloured by polymer | `x`, `R`, `pol_id` |
| [`LPlot.polymers_3d!`](@ref) | spheres, coloured by polymer (`plotdim = 3`) | `x`, `R`, `pol_id` |
| [`LPlot.ellipses!`](@ref) | ellipses from a shape tensor, rotated along the polarity | `x`, `p`, `Lambda` |
| [`LPlot.sphere!`](@ref) | a white background sphere, for particles confined to a sphere (`plotdim = 3`) | `x` |

## Rigid bodies

| Plot function | Shows | Needs particle fields |
|---|---|---|
| [`LPlot.shape_disks!`](@ref) | circles at the extent points, coloured by id | `xe`, `re`, `id` |
| [`LPlot.shape_disks_orientation!`](@ref) | circles at the extent points, coloured by polarity angle | `xe`, `re`, `p` |
| [`LPlot.shape_disks_type!`](@ref) | circles at the extent points, coloured by type | `xe`, `re`, `type` |
| [`LPlot.shape_points!`](@ref) | points at the extent points, coloured by id | `xe`, `re`, `id` |

## Vectors and labels

| Plot function | Shows | Needs particle fields |
|---|---|---|
| [`LPlot.directors!`](@ref) | arrows along the polarity | `x`, `p` |
| [`LPlot.nematic_directors!`](@ref) | double-headed arrows along `p` and `-p` | `x`, `p` |
| [`LPlot.velocity_vectors!`](@ref) | arrows along the velocity | `x`, `v` |
| [`LPlot.annotate_v!`](@ref) | the speed of every particle as a text label | `x`, `v` |

## Fields

These plot the first field of the field state.

| Plot function | Shows | Needs field fields |
|---|---|---|
| [`LPlot.field_magnitude!`](@ref) | heatmap of `C` (range 0–2) with the bin edges | `bin_centers`, `lbin`, `C` |
| [`LPlot.GPUfield_magnitude!`](@ref) | as `field_magnitude!`, for single-precision or GPU fields | `bin_centers`, `lbin`, `C` |
| [`LPlot.field_log_magnitude!`](@ref) | heatmap of `log10(C)` (range -4 to -2) | `bin_centers`, `C` |
| [`LPlot.field_magnitude_wgrid!`](@ref) | heatmap of `C` (range 0–1) with grid lines at the bin centres | `bin_centers`, `C` |
| [`LPlot.potential!`](@ref) | surface plot of `C` as an energy landscape (`plotdim = 3`) | `bin_centers`, `C` |

## Reference

```@docs
LPlot
```

### Particles as points

```@docs
LPlot.points!
LPlot.type_points!
LPlot.director_points!
LPlot.trajectories!
LPlot.Swarmalators!
```

### Particles at their true size

```@docs
LPlot.disks!
LPlot.disks_type!
LPlot.disks_orientation!
LPlot.disks_nematic_orientation!
LPlot.disks_v_orientation!
LPlot.disks_vp_phase_difference!
LPlot.disks_vx!
LPlot.disks_uw!
LPlot.transparant_disks!
LPlot.sized_points!
LPlot.type_sized_points!
LPlot.polymers!
LPlot.polymers_3d!
LPlot.ellipses!
LPlot.sphere!
```

### Rigid bodies

```@docs
LPlot.shape_disks!
LPlot.shape_disks_orientation!
LPlot.shape_disks_type!
LPlot.shape_points!
```

### Vectors and labels

```@docs
LPlot.directors!
LPlot.nematic_directors!
LPlot.velocity_vectors!
LPlot.annotate_v!
```

### Fields

```@docs
LPlot.field_magnitude!
LPlot.GPUfield_magnitude!
LPlot.field_log_magnitude!
LPlot.field_magnitude_wgrid!
LPlot.potential!
```

### Helpers

```@docs
LPlot.ellipse
LPlot.angle2range
```
