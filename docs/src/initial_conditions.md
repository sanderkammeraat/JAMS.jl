# Initial conditions

The `Initial` module has helper functions that generate positions for common initial
conditions and confining geometries. Each function returns two vectors, `xs` and `ys`, with
the x- and y-coordinates in the xy-plane, centred on the origin. Use them to build the
particles of the initial state:

```julia
xs, ys = Initial.random_in_disk(500, 20.0)

initial_state = ParticleState([
    Particles.Polar(id=i, type=1, x=[xs[i], ys[i], 0.0], p=normalize([randn(), randn(), 0.0]))
    for i in eachindex(xs)
])
```

## Walls of fixed particles

JAMS has no built-in walls. Instead, you can build a wall from particles that no DOF evolver
acts on, so they never move. [`Initial.box`](@ref) and [`Initial.ring`](@ref) generate the
positions for such walls. For example, active particles (type `1`) inside a square box of
wall particles (type `2`):

```julia
Lx, Ly = 50.0, 50.0
xs, ys = Initial.box(1.5, Lx, Ly)            # wall particles about 1.5 apart
Nb = length(xs)
Nint = 1000

inner = [Particles.Polar(id=i, type=1,
            x=[rand(Uniform(-Lx/2 + 4, Lx/2 - 4)), rand(Uniform(-Ly/2 + 4, Ly/2 - 4)), 0.0],
            p=normalize([randn(), randn(), 0.0])) for i in 1:Nint]
wall  = [Particles.Polar(id=Nint + i, type=2, x=[xs[i], ys[i], 0.0], p=[1.0, 0.0, 0.0])
            for i in 1:Nb]
initial_state = ParticleState(vcat(inner, wall))

forces = (Forces.self_propulsion(1, 0.01), Forces.repulsive_soft_disk([1, 2], [1.0 1.0; 1.0 1.0]))
dofevolvers = (DOFevolvers.overdamped_xvf(1), DOFevolvers.overdamped_pqT_xyc(1))  # type 2 stays fixed
```

Make the simulation box slightly larger than the wall (e.g. `sizes = (Lx + 4, Ly + 4, 1.0)`)
so the wall particles are inside it. Keep the wall particles closer together than the
diameter of the inner particles, so none can slip through. See `examples/confined_SABPs.jl`
(square box), `examples/confined_SABPs_disk.jl` (disk) and
`examples/confined_SABPs_annulus.jl` (annulus) for complete scripts.

## Available generators

| Function | Generates |
|---|---|
| [`Initial.box`](@ref) | positions along the edge of a rectangle, spaced about `l` apart |
| [`Initial.ring`](@ref) | positions on a circle, spaced an arc length of about `l` apart |
| [`Initial.random_in_disk`](@ref) | `N` random positions spread uniformly over a disk or annulus |
| [`Initial.stacked_polymers_at_angle`](@ref) | straight polymers filling a periodic box at a given packing fraction, with their polymer ids, positions along the polymer and the box size |

## Polymers

[`Initial.stacked_polymers_at_angle`](@ref) returns more than positions: also the polymer
id and position along the polymer of every monomer, and the box size `L` that gives the
requested packing fraction. Use them to build [`Particles.PolarPolymer`](@ref)s:

```julia
N_in_pol = 10
x, y, radii, pol_ids, ids_in_pol, L, Npols = Initial.stacked_polymers_at_angle(N_in_pol, 1.0, 0.9, 0.75, 80.0)

initial_state = ParticleState([
    Particles.PolarPolymer(id=i, type=1, pol_id=pol_ids[i], id_in_pol=ids_in_pol[i], pol_N=N_in_pol,
        R=radii[i], x=[x[i], y[i], 0.0], p=normalize([randn(), randn(), 0.0]))
    for i in 1:Npols*N_in_pol
])
sizes = (L, L, 2.0)   # with Periodic = true
```

See `examples/polymers.jl` for a complete script.

## Reference

```@docs
Initial.box
Initial.ring
Initial.random_in_disk
Initial.stacked_polymers_at_angle
```
