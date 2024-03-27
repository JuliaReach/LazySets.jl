# Optional Features

This section summarizes each optional dependency and the features available when used together with `LazySets.jl`.
Internally, optional dependencies in Julia are handled with the package [Requires.jl](https://github.com/JuliaPackaging/Requires.jl).

```@contents
Pages = ["optional_dependencies.md"]
```

## Installing all dependencies

Use the following command to install *all* optional dependencies. Installing all optional dependencies is required
if you want to run the full test suite and build the documentation locally.

```julia
julia> import Pkg; Pkg.add(["CDDLib",
                            "Distributions",
                            "Documenter",
                            "Expokit",
                            "ExponentialUtilities",
                            "IntervalMatrices",
                            "Makie",
                            "Optim",
                            "Polyhedra",
                            "RecipesBase",
                            "StaticArrays",
                            "Symbolics",
                            "TaylorModels"])
```

## Documentation

To build the project's documentation locally you need to install [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl).

|Dependency|Features|
|----------|-------|
|[`Documenter`](https://github.com/JuliaDocs/Documenter.jl) |Building the package's documentation.|
|[`LaTeXStrings`](https://github.com/stevengj/LaTeXStrings.jl) |Input and display of LaTeX equation strings, used in plots.|

## Exponential backends

Exponential backends are used to compute the action of matrix exponentials over sets lazily, e.g. ``\rho(d, e^{A \delta} X)`` for ``A`` large and sparse.

|Dependency|Features|
|----------|-------|
|[`ExponentialUtilities`](https://github.com/SciML/ExponentialUtilities.jl) | Utility functions for exponential integrators from the SciML scientific machine learning ecosystem.|
|[`Expokit`](https://github.com/acroy/Expokit.jl) |Julia implementation of EXPOKIT routines.|

## Interval methods

For validated numerics, we build upon the Julia ecosystem [JuliaIntervals](https://github.com/JuliaIntervals). The package [IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl), which is a hard dependency of LazySets, implements rigorous floating-point calculations using interval arithmetic in Julia and is the basis for the implementation of [`Interval`](@ref). There are other interval packages that can also be used in conjunction with LazySets and provide additional functionality.

The package [`IntervalLinearAlgebra`](https://github.com/JuliaIntervals/IntervalLinearAlgebra.jl) is not an optional dependency of LazySets, but these packages can be used in conjunction to compute [solution sets of interval linear systems](https://juliaintervals.github.io/IntervalLinearAlgebra.jl/dev/explanations/solution_set/).

|Dependency|Features|
|----------|-------|
|[`IntervalConstraintProgramming`](https://github.com/JuliaIntervals/IntervalConstraintProgramming.jl) |Conservative polyhedral approximation of a region defined implicitly via nonlinear constraints.|
|[`IntervalMatrices`](https://github.com/JuliaReach/IntervalMatrices.jl) | Set operations that involve matrices whose coefficients are intervals.|
|[`TaylorModels`](https://github.com/JuliaIntervals/TaylorModels.jl) |Taylor expansion of functions with rigorous interval remainder.|

## Optimization algorithms

Some computations require use of external numerical optimization solvers. The modeling language [`JuMP`](https://github.com/jump-dev/JuMP.jl) is loaded by default, together with the [GLPK](https://en.wikipedia.org/wiki/GNU_Linear_Programming_Kit) solver for linear programs (LPs). Other solvers can be loaded on-demand, even commercial ones, provided that you have the appropriate license. See JuMP's [documentation page on supported solvers](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers) for further details.

For other uses, such as line search methods to compute the support function of lazy intersections of certain sets, [`Optim`](https://github.com/JuliaNLSolvers/Optim.jl) can be optionally loaded.

|Dependency|Features|
|----------|-------|
|[`Optim`](https://github.com/JuliaNLSolvers/Optim.jl) |Optimization package in pure Julia.|


## Polyhedral computations

The optional package [Polyhedra.jl](http://github.com/JuliaPolyhedra/Polyhedra.jl) is required whenever you want to work with concrete polyhedra representations in dimension higher than 2, e.g. to solve the vertex enumeration problem (converting from constraint to vertex representation of a polytope). While Polyhedra implements its own default backend, it is also possible to load external ones, such as [cdd](https://www.inf.ethz.ch/personal/fukudak/cdd_home/) through [CDDLib.jl](https://github.com/JuliaPolyhedra/CDDLib.jl). For a list of available polyhedra backends, see the website of the [JuliaPolyhedra ecosystem](https://juliapolyhedra.github.io/).

|Dependency|Features|
|----------|-------|
|[`Polyhedra`](http://github.com/JuliaPolyhedra/Polyhedra.jl) |Concrete polyhedra library.|
|[`CDDLib`](https://github.com/JuliaPolyhedra/CDDLib.jl) |Polyhedral computations backend.|

## Random sampling and triangulation

|Dependency|Features|
|----------|-------|
|[`Distributions`](https://github.com/JuliaStats/Distributions.jl) |Random sampling.|
|[`MiniQhull`](https://github.com/gridap/MiniQhull.jl) |Voronoi-Delaunay triangulation of `LazySets` types.|


## Special array types

|Dependency|Features|
|----------|-------|
|[`StaticArrays`](https://github.com/JuliaArrays/StaticArrays.jl) |Statically defined arrays.|


## Symbolic utilities

|Dependency|Features|
|----------|-------|
|[`Symbolics`](https://github.com/JuliaSymbolics/Symbolics.jl) |Create sets using symbolic expressions.|

## Visualization

As in other Julia packages, the user has to manually install and load the package [`Plots`](https://github.com/JuliaPlots/Plots.jl) in order to produce visual results. As a plotting backend we recommend using [GR](https://github.com/jheinen/GR.jl), which is the default backend and is automatically installed together with Plots.

|Dependency|Features|
|----------|-------|
|[`Makie`](https://github.com/JuliaPlots/Makie.jl) |Mainly for 3D and interactive plots.|
|[`Plots`](https://github.com/JuliaPlots/Plots.jl) |Mainly 2D plots.|
