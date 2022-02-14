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
although some extra functionality is available provided that you load that package. 
For example, if you want to work with sets defined using simple algebraic expressions you can install
[Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) as usual with the package manager,
`] add Symbolics`, then load it together with `LazySets` to have new functionality.

```julia
julia> using LazySets, Symbolics

julia> var = @variables x[1:10]
(Num[x₁, x₂, x₃, x₄, x₅, x₆, x₇, x₈, x₉, x₁₀],)

julia> Hyperplane(x[1] + x[2] == 1/2, var)
Hyperplane{Float64,Array{Float64,1}}([1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 0.5)
```
defines the hyperplane $H : \{ x \in \mathbb{R}^{10} : x_1 + x_2 = 1/2\}$. The `Hyperplane` constructor
(and other constructors) automatically work with Symbolics's variables once `using Symbolics` is ran
in your session.

The full list of optional dependencies can be found in section [Optional Features](@ref).

## Running the unit tests

To run the full test suite, do

```julia
] test LazySets
```
Please note that it is required that you first install all the optional dependencies as specified in [Installing all dependencies](@ref).
It is possible to use optional flags to select a portion of tests to be run. See the script `test/runtests.jl` for details.
