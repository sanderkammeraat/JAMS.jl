# Forces

Forces add to the total force `f` and torque `T` of particles every timestep. They never
move particles themselves; that is what the  [DOF evolvers](dofevolvers.md) do.

Put the forces of a simulation in a tuple and pass it to the [`System`](@ref):

```julia
forces = (
    Forces.self_propulsion(1, 0.2),
    Forces.planar_rotational_noise(ontypes=1, Dr=0.01),
    Forces.repulsive_soft_disk(1, 1.0),
)
```

Most forces can be constructed with positional arguments in the order of their fields, and
those with `(; ...)` in their signature also with keywords.

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
| [`Forces.translational_noise`](@ref) | Gaussian white-noise on x, y and z. |

## Pair forces

| Force | Description |
|---|---|
| [`Forces.repulsive_soft_disk`](@ref) | Harmonic repulsion between overlapping particles, `U = k/2 (R_i + R_j - r)^2`. |
| [`Forces.morse`](@ref) | Morse potential: repulsive below contact, attractive beyond it, up to `rcut_pair_global`. |
| [`Forces.pairAN`](@ref) | Force (and optional torque) from the active nematic stresses of two nearby particles. |
| [`Forces.pair_polar_alignment`](@ref) | Torque that aligns the polarities of nearby particles (`p_i` towards `p_j`). |
| [`Forces.pair_nematic_alignment`](@ref) | Torque that aligns the polarities of nearby particles nematically (`p_i` towards `p_j` or `-p_j`). |

The aligners set their own range, `rcut` or `rfact * (R_i + R_j)`. Pairs are still
only evaluated within `rcut_pair_global`, so set it at least as large as that range.

### Polymer forces
The polymer forces have been developed as part of G. Martin's MSc thesis.
These pair forces need a particle type that knows which polymer it belongs to, such as
[`Particles.PolarPolymer`](@ref) (fields `pol_id`, `id_in_pol` and `pol_N`). A minimal
polymer model combines bonds, bending stiffness and excluded volume:

```julia
forces = (
    Forces.polymer_harmonic_stretch(1, false, 3.0, 0.75),  # bonds, length 0.75 * (R_i + R_j)
    Forces.polymer_harmonic_bend(1, 3.0),                  # bending stiffness
    Forces.polymer_repulsive_soft_disk(1, 1.0),            # soft disk repulsion, but not within bonds
)
dofevolvers = (DOFevolvers.overdamped_xvf(1), DOFevolvers.polymer_p_set(1))
```

[`DOFevolvers.polymer_p_set`](@ref) keeps the polarity `p` of each monomer along the
polymer, which the active polymer forces use. Bending acts between monomers two places
apart along a polymer, so set `rcut_pair_global` to more than twice the bond length. See
`examples/polymers.jl` for a complete script.

| Force | Description |
|---|---|
| [`Forces.polymer_harmonic_stretch`](@ref) | Harmonic bond between consecutive monomers of a polymer. |
| [`Forces.polymer_harmonic_bend`](@ref) | Bending forces. |
| [`Forces.polymer_repulsive_soft_disk`](@ref) | Soft-disk repulsion between all monomers except bonded ones. |
| [`Forces.polymer_pair_polar_nematic`](@ref) | Active sliding force between monomers of different polymers, set by their relative polarity. |
| [`Forces.polymer_pairAN`](@ref) | [`Forces.pairAN`](@ref) between monomers of different polymers only. |

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
Forces.translational_noise
```

### Pair forces

```@docs
Forces.repulsive_soft_disk
Forces.morse
Forces.pairAN
Forces.pair_polar_alignment
Forces.pair_nematic_alignment
```

### Polymer forces

```@docs
Forces.polymer_harmonic_stretch
Forces.polymer_harmonic_bend
Forces.polymer_repulsive_soft_disk
Forces.polymer_pair_polar_nematic
Forces.polymer_pairAN
```

### Abstract types

```@docs
Forces.Force
Forces.ExternalForce
Forces.PairForce
Forces.FieldForce
```
