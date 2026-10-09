# Extending JAMS

JAMS is meant to be easily extended. A new
force or DOF evolver consists of the combination of a struct(containing the parameters) plus one method (actually calculating the force). Thanks to Julia's multiple dispatch, the engine picks up the method.

## How forces see a particle: `p_i`

Forces and DOF evolvers receive a particle as `p_i`, a mutable view of particle i of the particle state.

Vectors are immutable `SVector`s, so you always assign a new value (`p_i.f += ...`) instead
of changing an element in place (e.g. `p_i.f[1] += ...` does not work).

## A custom external force

An external force acts on one particle at a time. Define a subtype of
[`Forces.ExternalForce`](@ref) and add a method to `Forces.contribute_external_force!`.

As an example, a constant force such as gravity needs to add a struct and a method to contributing external forces:

```julia
struct constant_force <: Forces.ExternalForce
    ontypes::Int
    F::SVector{3, Float64} #The constant force vector we are adding
end

function contribute_external_force!(p_i, t, dt, rngs_particles, system, force::constant_force)
    if p_i.type in force.ontypes
        p_i.f += force.F
    end
    return p_i
end
```

## A custom pair force

A pair force acts between two particles closer than `rcut_pair_global`. Define a subtype of
[`Forces.PairForce`](@ref) and add a method to `Forces.contribute_pair_force!`. The engine
passes both particles, the separation vector `dx = x_j - x_i` (the nearest periodic image
in a periodic box) and its length `dxn`.


The engine calls this method once for `(i, j)` and once for `(j, i)`, possibly on
different threads at the same time. So:

- only change `p_i`, never `p_j`
- compute the force *on* `i`: `dx` points from `i` to `j`, so a repulsive force has a
  minus sign

## A custom DOF evolver

A local DOF evolver updates one particle from its accumulated force and torque. Define a
subtype of [`DOFevolvers.LocalDOFevolver`](@ref) and add a method to
`DOFevolvers.evolve_locally!`. It must reset the forces or torques it used, otherwise they
keep adding up over the timesteps.

## A custom global DOF evolver

A global DOF evolver needs the state of other particles, for example its neighbours.
Define a subtype of [`DOFevolvers.GlobalDOFevolver`](@ref) and add a method to
`DOFevolvers.evolve_globally!`. It receives the current particle and field state and the Verlet
neighbour list, and returns the updated particle and field state.

## A custom particle type

If you need extra degrees of freedom or parameters, define your own particle type. It must
have at least these fields, which the engine uses:

| Field | Type | Used for |
|---|---|---|
| `id` | `Int64` | random number generator of the particle; ids must run from `1` to `N` |
| `type` | `Int64` | Applying forces and DOFevolvers only to certain types via `ontypes` of forces and DOF evolvers |
| `x` | `SVector{3,Float64}` | position, boundary conditions and neighbour search |
| `xuw` | `SVector{3,Float64}` | unwrapped position, set to `x` at `t = 0` |
| `f` | `SVector{3,Float64}` | total force, reset at `t = 0` |
| `T` | `SVector{3,Float64}` | total torque, reset at `t = 0` |

The polymer forces and [`DOFevolvers.polymer_p_set`](@ref) also need `pol_id`,
`id_in_pol` and `pol_N`, as in [`Particles.PolarPolymer`](@ref).

## A custom save function

A save function is called every `Tsave` timesteps with the HDF5 group of the current frame.
Write datasets to it by assigning arrays:

```julia
function save_positions!(current_frame_group, current_particle_state, current_field_state, n, Tsave, t, framecounter)
    current_frame_group["t"] = t
    current_frame_group["x"] = [p_i.x[1] for p_i in current_particle_state]
    return current_frame_group
end
```

Pass it as `save_functions = (save_positions!,)`. See [`Save.polar_particle!`](@ref) for
the built-in one.

## A custom plot function

A plot function is called once, when the live plot is set up. It receives the Makie figure
and axis, and `Observable`s of the particle and field state. Use `@lift` to derive the
data to plot, so that the plot updates whenever the state changes:

```julia

function points!(f, ax, cpsO, cfsO)
    x = @lift([p.x[1] for p in $cpsO])
    y = @lift([p.x[2] for p in $cpsO])
    scatter!(ax, x, y)
    return ax
end
```


## Tips for fast custom code

The force and evolver methods run for every particle (and every neighbour) in every
timestep, so small inefficiencies add up.

- **Use concrete field types in your structs**, or type parameters:
  `struct my_force{T} <: Forces.PairForce; k::T; end`. Avoid `Any`, `Real` or untyped
  fields.
- **Use `SVector` for vectors**, not `Vector`. `SVector` arithmetic does not allocate
  memory; `Vector` arithmetic allocates on every call.
- **Put forces and evolvers in tuples**, not vectors, as in all examples. The engine then
  compiles one specialized loop for your exact combination.
- **Do not use global variables** inside forces; store parameters in the struct instead.
- **Keep `rcut_pair_global` as small as your pair forces allow** (see
  [Pair cutoff and neighbour search](@ref)).
