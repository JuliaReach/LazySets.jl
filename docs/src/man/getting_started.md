# Getting Started

In this section we review the recommended setup to start working with this
package.

```@contents
Pages = ["getting_started.md"]
```

## Setup

This package requires Julia v1.0 or later.
Refer to the [official documentation](https://julialang.org/downloads/) on how
to install it for your system.
Below we explain the steps for setting up `LazySets` on your system and checking
that it builds correctly.


### Installation

To install `LazySets`, use the following command inside Julia's REPL:

```julia
julia> import Pkg; Pkg.add("LazySets")
```
or replace `add` by `clone` if you want to develop the code.
The full list of dependencies (which are automatically installed) is specified
in the [Project.toml](https://github.com/JuliaReach/LazySets.jl/blob/master/Project.toml) file.


### Building the package

Use the following command from Julia's REPL:
```
julia> using LazySets
```
This should precompile the package and make it available afterward.


## Workflow tips

There are different ways to use Julia: from the terminal, from the Julia REPL,
from IJulia (i.e., Jupyter notebook), from Juno, etc.
If you do not have a preferred choice, we recommend using `LazySets` through
IJulia;
one reason is that the visualization is conveniently embedded into the notebook,
and it can be exported into different formats, among other benefits.
On the other hand, for development purposes you probably want to use the REPL or
the Juno environment.


## Updating

After working with `LazySets` for some time, you may want to get the newest
version.
For this you can use the following command (e.g., from the REPL):
```julia
Pkg.checkout("LazySets")
```
That will check out the latest version in the `master` branch, and precompile it
the next time you enter a session and execute `using LazySets`.

## Optional dependencies

An optional dependency is a package that is not required to compile and use `LazySets.jl`,
although some extra functionality is available provided that you load that package. Internally, optional
dependencies in Julia are handled with the package [Requires.jl](https://github.com/JuliaPackaging/Requires.jl).

For example, if you want to work with sets defined symbolically as in `var = @variables x[1:10]; Hyperplane(x[1] + x[2] == 1/2, var)` to define the hyperplane $H : \{ x \in \mathbb{R}^{10} : x_1 + x_2 == 1/2\}$, install [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) so that such constructor (and other symbolic features) are automatically available once doing `using ModelingToolkit` in your session.
Another optional dependency is [`Expokit.jl`](https://github.com/acroy/Expokit.jl), a package that provides
lazy matrix exponentiation routines for linear maps of `LazySets` types.

The following table summarizes each optional dependency and the features available together with `LazySets.jl`:

|Dependency|Features|
|----------|-------|
|`CDDLib` |Polyhedral computations backend|
|`Distributions` |Random sampling|
|`Documenter` |Building the projects' documentation|
|`Expokit` |Computing the action of a matrix exponential over a set|
|`GR` |Plotting backend|
|`IntervalConstraintProgramming` |Conservative polyhedral approximation of a region defined implicitly via nonlinear constraints|
|`IntervalMatrices` |Set operations that involve matrices whose coefficients are intervals|
|`Makie` |Plotting library, mainly for 3D and interactive plots|
|`MiniQhull` |Voronoi-Delaunay triangulation of `LazySets` types|
|`ModelingToolkit` |Create sets using symbolic expressions|
|`Optim` |Numerical optimization package|
|`Plots` |Plotting library, mainly 2D plots|
|`Polyhedra` |Concrete polyhedra library|
|`StaticArrays` |Statically defined arrays, mainly for low dimensions|
|`TaylorModels` |Taylor expansion of functions with rigorous interval remainder|

Use the following command to install *all* optional dependencies:

```julia
julia> ] add CDDLib Distributions Documenter Expokit GR IntervalConstraintProgramming IntervalMatrices Makie MiniQhull ModelingToolkit Optim Plots Polyhedra StaticArrays TaylorModels
```
