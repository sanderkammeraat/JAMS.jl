# Getting started

## Installation

In a Julia script or in the Julia REPL, run

```julia
using Pkg
Pkg.add(url="https://github.com/sanderkammeraat/JAMS.jl")
```

and load the package with

```julia
using JAMS
```

`using JAMS` also loads `Random`, `Distributions` and `LinearAlgebra`, so functions such as
`rand(Uniform(0, 1))` and `normalize` are available straight away.

JAMS runs in parallel on as many threads as Julia is started with. If one uses VS Code, the number of threads used by Julia are set automatically.  Running from the terminal, you have to specify the number of threads manually. To start Julia with, for example, 4 threads, do:

```bash
julia -t 4
```

## An example simulation: Active Brownian particles

This example simulates 1000 active Brownian particles (ABPs) in a periodic 2D box. Each
particle moves with a constant speed along its polarity vector, the polarity diffuses
randomly, and overlapping particles repel each other according to a soft disk force.

### 1. Choose the forces

```julia
using JAMS 
forces = (
    Forces.self_propulsion(1, 0.2),                      # speed v0 = 0.2
    Forces.planar_rotational_noise(ontypes=1, Dr=0.01),  # rotational diffusion
    Forces.repulsive_soft_disk(1, 1.0),                  # harmonic repulsion, stiffness 1
)
```

The first argument of every force, `ontypes`, lists the particle types it acts on. Here all
particles have type `1`. As you can see, the definition for forces in JAMS is a bit loose: for example the rotational 
noise is strictly a torque instead of a force.

### 2. Choose the DOF evolvers

A DOF (degree of freedom) evolver takes the total force and torque and applies it to a particle. For overdamped particles in the xy-plane, use one evolver for the position and one
for the polarity:

```julia
dofevolvers = (
    DOFevolvers.overdamped_xvf(1),       # x += f/zeta * dt
    DOFevolvers.overdamped_pqT_xyc(1),   # rotate p by T/zeta_R * dt about z
)
```

### 3. Setup the initial particle condition

```julia
N = 1000
ϕ = 1.0                                   # packing fraction
poly = 0.15                               # polydispersity of the radii
Rs = rand(Uniform(1 - poly, 1 + poly), N)
L = sqrt(pi * sum(Rs .^ 2) / ϕ)           # box size for packing fraction ϕ

initial_state = ParticleState([
    Particles.Polar(
        id = i,
        type = 1,
        R = Rs[i],
        x = [rand(Uniform(-L/2, L/2)), rand(Uniform(-L/2, L/2)), 0.0],
        p = normalize([randn(), randn(), 0.0]),
    )
    for i in 1:N
])
```

Positions always have three components. For a 2D simulation, set the z-component to `0`.
Particle ids must run from `1` to `N`.

### 4. Define the system

```julia
system = System(
    sizes = (L, L, 2.0),
    initial_particle_state = initial_state,
    forces = forces,
    dofevolvers = dofevolvers,
    Periodic = true,
    rcut_pair_global = 2 * (1 + poly),     # largest possible R_i + R_j
)
```

The box runs from `-L/2` to `L/2` in each direction. `rcut_pair_global` is the largest
distance at which pair forces are evaluated. Soft disks only touch below `R_i + R_j`, so
twice the largest radius is enough; a larger cutoff gives the same result but is slower.

### 5. Run it

```julia
sim = Euler_integrator(system, 0.01, 100.0; seed=42)
```

This integrates from `t = 0` to `t = 100` with timestep `dt = 0.01`, showing a progress bar.
The result is a [`JAMS.SIM`](@ref):

```julia
sim.final_particle_state.x        # final positions of all particles
sim.final_particle_state[1]       # final state of particle 1
```

With the same `seed`, the simulation gives the same result every time, also with a
different number of threads.

## Saving data

To write data to disk, give a save interval, one or more save functions and a folder:

```julia
sim = Euler_integrator(system, 0.01, 100.0;
    seed = 42,
    Tsave = 100,                              # save every 100 timesteps
    save_functions = (Save.polar_particle!,),
    save_folder_path = joinpath(pwd(), "ABP_example"),
)
```

This creates `ABP_example/raw_data.h5` with one group per saved frame (`frames/1`,
`frames/2`, ...). Read it back with HDF5.jl or more conveniently in Julia using JLD2:

JAMS prevents writing into a folder that already contains output files, so an earlier run
is never overwritten. Use a new folder, or set `save_tag` to give the files a prefix. See
[`Euler_integrator`](@ref) for all output files.

## Live plotting

For live plotting, load GLMakie *after* JAMS:

```julia
using JAMS
using GLMakie
```
This enables the optional live plotting extension of JAMS and compiles it.
and pass a plot interval and plot functions:

```julia
sim = Euler_integrator(system, 0.01, 1000.0;
    Tplot = 10,                                                #update plot every 10 time steps
    plot_functions = (LPlot.disks_orientation!, LPlot.directors!),
)
```

This opens a window that shows the particles as disks coloured by their orientation, with
arrows along their polarity. Closing the window stops the simulation. To record a video,
also set `record_folder_path`. The [`LPlot`](@ref) module lists all plot functions.

On a cluster without a display, leave out `using GLMakie` and `Tplot`.

## Next steps

- [Overview](overview.md) explains the building blocks and the order of operations in a timestep.
- [Extending JAMS](@ref) shows how to write your own forces and update rules.
- The `examples/` folder has complete scripts, including the one above with self-alignment
  (`examples/ABPs.jl`).
