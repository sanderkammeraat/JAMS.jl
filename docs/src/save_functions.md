# Save functions

Save functions write the particle (and field) state to disk during a simulation. Pass a
tuple of them to [`Euler_integrator`](@ref), together with a save interval and a folder:

```julia
sim = Euler_integrator(system, 0.01, 100.0;
    Tsave = 100,                                  # save every 100 timesteps
    save_functions = (Save.polar_particle!,),
    save_folder_path = joinpath(pwd(), "run1"),
    save_tag = "abp",                             # optional prefix for the file names
)
```

Every `Tsave` timesteps, JAMS creates a new group `frames/1`, `frames/2`, ... in the HDF5
file `raw_data.h5` and calls each save function with that group. The save function decides
which datasets to write. The saved positions are those at the start of the timestep, before
they are updated.

Besides `raw_data.h5`, JAMS writes `JAMs_container.jld2` (the [`System`](@ref) and
integration settings) and `JAMs_final_state.jld2` (the final state). With `save_tag`, all
three file names get the tag as a prefix, e.g. `abp_raw_data.h5`. JAMS refuses to start if
these files already exist in the folder, so earlier results are never overwritten.

## Reading the data back

```julia
using HDF5

h5open(joinpath("run1", "abp_raw_data.h5"), "r") do file
    nframes = length(file["frames"])
    for k in 1:nframes
        frame = file["frames/$k"]
        t = read(frame["t"])
        x = read(frame["x"])
        y = read(frame["y"])
        # ... analyse frame k
    end
end
```

The file also contains the parameters of the system, forces and DOF evolvers, and the
integration settings (under `system` and `integration_info`), so each output file describes
the simulation that produced it.

## Available save functions

| Save function | Saves | For particles with fields |
|---|---|---|
| [`Save.polar_particle!`](@ref) | `n`, `t`, `id`, `type`, `R`, position, unwrapped position, velocity, polarity and angular velocity, split into x-, y- and z-components | `id`, `type`, `R`, `x`, `xuw`, `v`, `p`, `q` |

To save other quantities, or for your own particle type, see
[A custom save function](@ref).

## Reference

```@docs
Save.polar_particle!
```
