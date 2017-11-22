# Getting Started

In this section we review the recommended setup to start working with this package.

```@contents
Pages = ["getting_started.md"]
```

## Setup

This package requires Julia v0.6 or later. Refer to the official documentation
on how to install it for your system. Below we explains the steps for setting up
`LazySets` in your system and checking that it builds correctly.

### Installation

To install `LazySets` use the following command inside Julia's REPL:

```julia
Pkg.clone("https://github.com/JuliaReach/LazySets.jl")
```
The dependencies of `LazySets`, such as [Expokit](https://github.com/acroy/Expokit.jl)
-- which provides lazy matrix exponentiation routines -- are automatically installed
through Julia's package manager. The full list of dependencies is specified in the
`REQUIRE` file. The minimal Julia version required is v0.6.0. 

### Testing

Unit tests execute specific portions of the library code, checking that the result
produced is the expected one. To run the unit tests run the following command
in the terminal:

```julia
$ julia --color=yes test/runtests.jl
```

Alternatively, you can test the package in Julia's REPL with the command:

```julia
julia> Pkg.test("LazySets")
```

## Workflow tips

There are different ways to use Julia: from the terminal (so called REPL), from
IJulia (i.e. Jupyter notebook), from Juno, ... If you don't have a preferred
choice, we recommend using `LazySets` through IJulia; one reason is that the
visualization is conveniently embedded into the notebook, and it can be exported
into different formats, among other benefits. On the other hand, for development
purposes you'll probably prefer using the REPL or the Juno environment.

## Updating

After working with `LazySets` for some time, you may want to get the newest version.
For this you can use this command:
```julia
Pkg.checkout("LazySets")
```
That will checkout the latest version `master` branch, and precompile it the next
time you enter a session and do `using LazySets`.
