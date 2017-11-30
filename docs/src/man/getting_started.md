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
Pkg.clone("https://github.com/JuliaReach/LazySets.jl")
```
The dependencies of `LazySets`, such as
[Expokit.jl](https://github.com/acroy/Expokit.jl)
-- which provides lazy matrix exponentiation routines -- are automatically
installed through Julia's package manager.
The full list of dependencies is specified in the `REQUIRE` file.


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


## Development

If you want to contribute to this package, you are very welcome.

To preserve maintainability of this package, we have three policies:
pull requests, documentation, and unit testing.


### Pull requests

We use a standard pull request policy: You work in a private branch and
eventually add a pull request, which is then reviewed by other programmers and
merged into the `master` branch.


### Documentation

New functions and types should be documented according to our
[guidelines](https://github.com/JuliaReach/LazySets.jl/wiki/Documentation-Guidelines)
directly in the source code.

You can view the source code documentation from inside the REPL by typing `?`
followed by the name of the type or function.
For example, the following command will print the documentation of the `LazySet`
type:

```julia
julia> ?LazySet
```

To create this documentation you are currently reading, we use
[Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/).
The sources for creating this documentation are found in `docs/src`.
You can easily include the documentation that you wrote for your functions or
types there (see the
[Documenter.jl guide](https://juliadocs.github.io/Documenter.jl/stable/man/guide/)
or our sources for examples).

To generate the documentation locally, execute the following command in the
terminal:

```
$ julia --color=yes docs/make.jl
```


### Unit testing

Unit tests execute specific portions of the library code and check that the
produced result is the expected one.

When you modify code in this package, you should make sure that all unit tests
pass.
We also advise adding new unit tests when adding new features to ensure
long-term support of your contributions.

Thanks to continuous integration machinery, all unit tests are automatically
executed for every branch (and hence also every pull request) on the server.

To run the unit tests locally, execute the following command in the terminal:

```
$ julia --color=yes test/runtests.jl
```

Alternatively, you can achieve the same from inside the REPL using the following
command:

```julia
julia> Pkg.test("LazySets")
```
