# JAMS.jl

*Jamming and Active Matter Simulations in Julia*

JAMS is a Julia package to simulate (Soft) Active Matter. Its design is modular, so that you can mix and match different forces, particles and fields.

The purpose of this package is twofold: on one hand it provides a convenient way to explore new active matter models by providing flexible construction of (types of) forces and particles (e.g. simple polar particles or polymers). The exploration is facilated by an optional live plotting extension, leveraging GLMakie's efficient GPU plotting to render e.g. particle's positions, velocity vectors or polarties to quickly gauge what the system behaves like for different parameter values.

The second is to be performant to run production simulations for actual scientific analysis. The package has been through extensive profiling, is multi-threaded and easily runs on head-less clusters. Output is stored in the HDF5 format.

## Where to start

- [Getting started](@ref) installs JAMS and walks through a first simulation of Active
  Brownian particles.
- [Overview](overview.md) provieds an overview of the JAMS code, specifically how particles, forces, DOF evolvers and the `System` struct fit together.
- [Particles](particles.md), [Forces](forces.md), [DOF evolvers](dofevolvers.md), [Live plotting](live_plotting.md), [Save functions](save_functions.md) and [Initial conditions](initial_conditions.md) list all available particle types, forces, DOF evolvers, live plotting functions, save functions and initial condition generators.
- [Extending JAMS](@ref) shows how new particles, forces and DOFevolvers etc. can be added.
- [API reference](@ref) lists every type and function.

The `examples/` folder of the repository contains complete scripts, such as active
Brownian particles, active particles in a soft gel, and confined self-aligning particles.

## Table of contents

```@contents
Pages = ["getting_started.md", "overview.md", "particles.md", "forces.md", "dofevolvers.md", "live_plotting.md", "save_functions.md", "initial_conditions.md", "extending.md", "api.md"]
Depth = 2
```
