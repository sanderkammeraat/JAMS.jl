# DOF evolvers

A DOF (degree of freedom) evolver turns the force `f` and torque `T` accumulated during a
timestep into motion, for example `v = f/zeta` and `x += v * dt`, and then resets `f` and
`T` to zero for the next timestep.

Put the DOF evolvers of a simulation in a tuple and pass it to the [`System`](@ref):

```julia
dofevolvers = (DOFevolvers.overdamped_xvf(1), DOFevolvers.overdamped_pqT_xyc(1))
```

Like forces, every DOF evolver has an `ontypes` field with the particle type or types it
acts on. Particles that no DOF evolver acts on stay where they are, which is useful for
fixed walls.

Position and orientation are updated by separate evolvers, so a simulation usually has one
of each:

| Simulation | Position | Orientation |
|---|---|---|
| Overdamped particles in the xy-plane | [`DOFevolvers.overdamped_xvf`](@ref) | [`DOFevolvers.overdamped_pqT_xyc`](@ref) |
| Overdamped particles in 3D | [`DOFevolvers.overdamped_xvf`](@ref) | [`DOFevolvers.overdamped_pqT`](@ref) |
| Overdamped polymers, polarity along the polymer | [`DOFevolvers.overdamped_xvf`](@ref) | [`DOFevolvers.polymer_p_set`](@ref) |

To write your own, see [A custom DOF evolver](@ref).

## Local DOF evolvers

Local DOF evolvers update one particle at a time, using only that particle's own fields.

| DOF evolver | Description |
|---|---|
| [`DOFevolvers.overdamped_xvf`](@ref) | Overdamped position update: `v = f/zeta`, `x += v dt`, `xuw += v dt`, then resets `f`. |
| [`DOFevolvers.overdamped_pqT`](@ref) | Overdamped polarity update in 3D: rotates `p` with angular velocity `q = T/zeta_R`, then resets `T`. |
| [`DOFevolvers.overdamped_pqT_xyc`](@ref) | Overdamped polarity update in the xy-plane: exact rotation about z by `T[3]/zeta_R * dt`, then resets `T`. |

## Global DOF evolvers

Global DOF evolvers need the state of other particles as well, for example the positions
of neighbouring monomers in a polymer. They receive the whole particle state and the
neighbour list, and are by construction guaranteed to run after all local DOF evolvers.

| DOF evolver | Description |
|---|---|
| [`DOFevolvers.polymer_p_set`](@ref) | Sets the polarity `p` of each monomer along the local tangent of its polymer. |

## Field DOF evolvers

No field DOF evolvers are included yet.

## Abstract types

| Type | Description |
|---|---|
| [`DOFevolvers.DOFevolver`](@ref) | Supertype of all DOF evolvers. |
| [`DOFevolvers.LocalDOFevolver`](@ref) | Supertype of local DOF evolvers; describes the method a new one needs. |
| [`DOFevolvers.GlobalDOFevolver`](@ref) | Supertype of global DOF evolvers. |
| [`DOFevolvers.FieldDOFevolver`](@ref) | Supertype of field DOF evolvers. |

## Reference

```@docs
DOFevolvers.overdamped_xvf
DOFevolvers.overdamped_pqT
DOFevolvers.overdamped_pqT_xyc
DOFevolvers.polymer_p_set
DOFevolvers.DOFevolver
DOFevolvers.LocalDOFevolver
DOFevolvers.GlobalDOFevolver
DOFevolvers.FieldDOFevolver
```
