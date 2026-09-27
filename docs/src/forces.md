# Forces

Forces add to the total force `f` and torque `T` of particles every timestep. They never
move particles themselves; that is the job of the [DOF evolvers](dofevolvers.md).

Put the forces of a simulation in a tuple and pass it to the [`System`](@ref):

```julia
forces = (
    Forces.self_propulsion(1, 0.2),
    Forces.planar_rotational_noise(ontypes=1, Dr=0.01),
    Forces.repulsive_soft_disk(1, 1.0),
)
```

The first field of every force, `ontypes`, is the particle type (an `Int`) or types (a
`Vector{Int}`) it acts on. A pair force only acts between two particles if both have a type
in `ontypes`.

There are three kinds of forces:

- **External forces** act on one particle at a time.
- **Pair forces** act between two particles closer than `rcut_pair_global`.
- **Field forces** act between particles and fields.

To write your own, see [A custom external force](@ref) and [A custom pair force](@ref).

## External forces

| Force | Description |
|---|---|
| [`Forces.self_propulsion`](@ref) | Constant force `zeta * v0 * p` along the polarity, giving self-propulsion speed `v0`. |
| [`Forces.planar_rotational_noise`](@ref) | Gaussian white-noise torque about a fixed axis, giving rotational diffusion of the polarity. |
| [`Forces.self_align_with_v`](@ref) | Torque that rotates the polarity towards the particle's own velocity. |

## Pair forces

| Force | Description |
|---|---|
| [`Forces.repulsive_soft_disk`](@ref) | Harmonic repulsion between overlapping particles, `U = k/2 (R_i + R_j - r)^2`. |
| [`Forces.morse`](@ref) | Morse potential: repulsive below contact, attractive beyond it, up to `rcut_pair_global`. |

## Field forces

No field forces are included yet.

## Abstract types

| Type | Description |
|---|---|
| [`Forces.Force`](@ref) | Supertype of all forces. |
| [`Forces.ExternalForce`](@ref) | Supertype of external forces; describes the method a new external force needs. |
| [`Forces.PairForce`](@ref) | Supertype of pair forces; describes the method a new pair force needs. |
| [`Forces.FieldForce`](@ref) | Supertype of field forces. |

## Reference

### External forces

```@docs
Forces.self_propulsion
Forces.planar_rotational_noise
Forces.self_align_with_v
```

### Pair forces

```@docs
Forces.repulsive_soft_disk
Forces.morse
```

### Abstract types

```@docs
Forces.Force
Forces.ExternalForce
Forces.PairForce
Forces.FieldForce
```
