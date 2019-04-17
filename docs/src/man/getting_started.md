# Getting Started

In this section we review the recommended setup to start working with this
package.

```@contents
Pages = ["getting_started.md"]
```

## Setup

This package requires Julia v0.6 or later.
Refer to the [official documentation](https://julialang.org/downloads/) on how
to install it for your system.
Below we explain the steps for setting up `LazySets` on your system and checking
that it builds correctly.


### Installation

To install `LazySets`, use the following command inside Julia's REPL:

```julia
Pkg.add("LazySets")
```
or replace `add` by `clone` if you want to develop the code.
The full list of dependencies (which are automatically installed) is specified
in the `REQUIRE` file.


### Building the package

Use the following command from Julia's REPL:
```
julia> using LazySets
```
This should precompile the package and make it available afterward.


## Optional dependencies

Some optional dependencies, such as
[`Expokit.jl`](https://github.com/acroy/Expokit.jl) (a package that provides
lazy matrix exponentiation routines), are not installed by default.
When loading the corresponding packages in addition, new functionality in
`LazySets` is added automatically (a feature that is possible through the
[`Requires`](https://github.com/MikeInnes/Requires.jl) package).


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
