# Design

JAMS is modular thanks to its abstraction of the general building blocks that are naturally part of a simulation. This page explains
what each abstract building block is and how together they define a complete simulation in JAMS.

## Building blocks

| Building block | What it is | Examples |
|---|---|---|
| Particle type | An immutable struct with the fields of one particle | [`Particles.Polar`](@ref) |
| Particle state | All particles, stored column by column | [`ParticleState`](@ref) |
| Force | Adds to the force `f` and/or torque `T` of particles | [`Forces.self_propulsion`](@ref), [`Forces.repulsive_soft_disk`](@ref) |
| DOF evolver | Takes the total force and torque `f` and `T` of particle and evolves the particle, then resets them | [`DOFevolvers.overdamped_xvf`](@ref) |
| System | Defines the system: sets the system sizes, the initial state, the forces and the DOF evolvers | [`System`](@ref) |
| Integrator | Taking in System, it runs the simulation and returns the final state | [`Euler_integrator`](@ref) |

### Particles and the particle state

A particle type such as [`Particles.Polar`](@ref) is an immutable struct: `id`, `type`,
radius `R`, position `x`, polarity `p`, force `f`, and so on. Vectors are 3-component
`SVector`s from StaticArrays.jl, also in 2D simulations.

All particles together form a [`ParticleState`](@ref), which is a `StructArray`. It stores
each field as its own column: `state.x` is the vector of all positions, `state.R` the
vector of all radii. This is used for peformance.

### Particle types and `ontypes`

Every particle has an integer `type`. Every force and DOF evolver has an `ontypes` field
with the type (an `Int`) or types (a `Vector{Int}`) it acts on. This is how you give
different particles different physics. For example, in `examples/confined_SABPs.jl`
the active particles have type `1` and the pinned particles type `2`:

```julia
forces = (
    Forces.self_propulsion(1, 0.01),                   # only type 1 is active
    Forces.repulsive_soft_disk([1, 2], [1 2; 2 2]),    # all particles repel each other, k_11 =1, k_12 =2, k_21 =2, k_22=2 
)
dofevolvers = (DOFevolvers.overdamped_xvf(1),)         # only type 1 moves
```

Particles that no DOF evolver acts on never move, so the type-2 particles form a fixed boundary.

For pair forces, both particles must have a type in `ontypes`. Parameters that depend on
the pair of types, such as the stiffness of [`Forces.repulsive_soft_disk`](@ref), can be
given as a matrix indexed by the two types.

### Forces

There are currently three kinds of forces, each with its own abstract type:

- [`Forces.ExternalForce`](@ref): acts on one particle at a time, e.g. self-propulsion,
  noise or an external field. It can be calculated for each particle independently.
- [`Forces.PairForce`](@ref): acts between two particles closer than `rcut_pair_global`,
  e.g. soft repulsion. It needs the relative positions of two particles.
- [`Forces.FieldForce`](@ref): acts between particles and fields.

Forces only *add* to `f` and `T`. They never move particles themselves.

### DOF evolvers

A DOF (Degree Of Freedom) evolver updates a particle's state from the accumulated `f` and
`T`, for example `v = f/ζ` and `x += v dt` for an overdamped particle, and then sets `f` and
`T` back to zero for the next timestep. Because the forces and the equations of motion are
separate, you can reuse the same forces with overdamped or inertial dynamics.

A position evolver and an orientation evolver are independent, so a typical simulation
uses both, e.g. [`DOFevolvers.overdamped_xvf`](@ref) and
[`DOFevolvers.overdamped_pqT_xyc`](@ref).

## The box and boundary conditions

The box is set by `sizes = (Lx, Ly, Lz)`, must all be Floats (!). The system spans  `-L/2` to `L/2` in each direction.
Every particle must start inside it.

- **Periodic** (`Periodic = true`): all three directions are periodic, and pair forces use
  the nearest periodic image of each neighbour. Keep the box length in every direction
  the particles move in at least twice `rcut_pair_global`.
- **Closed** (`Periodic = false`): If a particle ends up outside the
  box, the simulation stops with an error.

**2D simulations** use the same 3D code: put all particles at `z = 0`, give the box a small
`Lz=1.` (it only needs to be positive), and use forces and evolvers that keep particles in the
plane, such as [`DOFevolvers.overdamped_pqT_xyc`](@ref) for the orientation.

## Pair cutoff and neighbour search

Pair forces are only evaluated for pairs closer than `rcut_pair_global`. To find those
pairs quickly, JAMS divides the box into cells at least `rcut_pair_global` wide, and only
compares each particle with particles in its own and the neighbouring cells. The cost per
timestep then grows linearly with the number of particles.

Choose `rcut_pair_global` as the range of your longest pair force, and no larger:

- for [`Forces.repulsive_soft_disk`](@ref): twice the largest radius
- for [`Forces.morse`](@ref): the distance where the attraction can be neglected

A cutoff that is too small silently cuts off interactions; one that is too large gives the
same result but costs more time.

## System: how everything comes together.

The `System` struct [] defines the complete system to be simulated, by collecting the system sizes, initial state, the forces, DOF evolvers and cutoff range.

## Euler_integrator: evolving the system over time

[`Euler_integrator`](@ref) loops over `t = 0, dt, 2dt, …, t_stop` and evolves the `System` struct over time. It specifies properties of the simulation during runtime, e.g. the plotting and saving. It returns the final state to easily allow chaining simulations.


## Chaining simulations

[`Euler_integrator`](@ref) returns a [`JAMS.SIM`](@ref) with the final state. You can use
it as the initial state of a next simulation, for example to relax a packing before
switching on activity:

```julia
relaxed = Euler_integrator(passive_system, 0.05, 1000.0)

active_system = System(
    sizes = relaxed.system.sizes,
    initial_particle_state = relaxed.final_particle_state,
    forces = active_forces,
    dofevolvers = dofevolvers,
    Periodic = true,
    rcut_pair_global = relaxed.system.rcut_pair_global,
)
active = Euler_integrator(active_system, 0.01, 1000.0)
```

Each simulation starts again at `t = 0`, so `f` and `T` are reset and the unwrapped
positions `xuw` restart from the current `x`.

## Units

JAMS has no built-in units: all quantities are numbers in whatever units you choose.

## Reproducibility and threads

JAMS runs the force, DOF evolver and boundary-condition loops in parallel, on the number of
threads Julia was started with e.g. (`julia -t 4`).

Every particle gets its own random number generator, seeded from the master seed and the
particle id. The random numbers a particle draws therefore do not depend on which thread
handles it, and a simulation with the same `seed` gives the same result on any number of
threads.

- Pass `seed` to [`Euler_integrator`](@ref) for a reproducible run. Without it, JAMS picks
  a random seed and stores it in the output files (`integration_info/master_seed`), so a
  saved run can still be repeated.
- `Euler_integrator` reseeds Julia's global random number generator with this seed. Create
  your initial conditions before calling it, or seed them separately with
  `Random.seed!`.
- In your own forces, draw random numbers only from `rngs_particles[p_i.id]`, never with a
  plain `rand()`, to keep this guarantee.