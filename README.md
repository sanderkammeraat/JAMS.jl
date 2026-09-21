
# JAMS: Jamming and Active Matter Simulations in Julia
[![Build Status](https://github.com/sanderkammeraat/JAMS.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/sanderkammeraat/JAMS.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Documentation Status](https://github.com/sanderkammeraat/JAMS.jl/actions/workflows/documentation.yml/badge.svg)](https://github.com/sanderkammeraat/JAMS.jl/actions/workflows/documentation.yml)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://sanderkammeraat.github.io/JAMS.jl/dev/)

<p align="center">
<img width="435" height="161" alt="JAMS_logo" src="https://github.com/user-attachments/assets/63f22d94-1de9-4556-a768-fd16407add82" />
</p>



## Description
This is a Julia package to simulate (Soft) Active Matter. Its design is modular, so that you can mix and match different forces, particles and fields.

The purpose of this package is twofold: on one hand it provides a convenient way to explore new active matter models by providing flexible construction of (types of) forces and particles (e.g. simple polar particles or polymers). The exploration is facilated by an optional live plotting extension, leveraging GLMakie's efficient GPU plotting to render e.g. particle's positions, velocity vectors or polarties to quickly gauge what the system behaves like for different parameter values.

The second is to be performant to run production simulations for actual scientific analysis. The package has been through extensive profiling, is multi-threaded and easily runs on head-less clusters. Output is stored in the HDF5 format.


## Installation
In a Julia script or from the REPL run
```
using Pkg
Pkg.add(url="https://github.com/sanderkammeraat/JAMS.jl")
```

You can then use it in a Julia script or REPL by importing the package:

```using JAMS```

## Optional Live plotting
Finally, to make use of live plotting, add the import of GLMakie *after* the import of JAMS, this will automatically load the live plotting extension:
```
using JAMS
using GLMakie
```
This will precompile the optional live plotting extension and enables it.

### Example
See ```examples/ABPs.jl``` to see an example of how the live plotting looks like. It is a simulation of Active Brownian particles with self-alignment, so that the systems starts to flock in periodic boundary conditions.
Here is the equivalent GIF showing a small system of 1000 particles to illustrate how the live plotting looks like.  

<img width="320" height="320" alt="output2" src="https://github.com/user-attachments/assets/af768ad2-6279-4b2e-9091-56187d0b38a5" />
