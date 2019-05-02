# About

This page contains some general information about this project, and
recommendations about contributing.

```@contents
Pages = ["about.md"]
```

## Contributing

If you like this package, consider contributing!
You can send bug reports (or fix them and send your code), add examples to the
documentation, or propose new features.

Below some conventions that we follow when contributing to this package are
detailed.
For specific guidelines on documentation, see the
[Documentations Guidelines wiki](https://github.com/JuliaReach/LazySets.jl/wiki/Documentation-Guidelines).

### Branches and pull requests (PR)

We use a standard pull request policy:
You work in a private branch and eventually add a pull request, which is then
reviewed by other programmers and merged into the `master` branch.

Each pull request should be pushed in a new branch with the name of the author
followed by a descriptive name, e.g., `mforets/my_feature`.
If the branch is associated to a previous discussion in one issue, we use the
name of the issue for easier lookup, e.g., `mforets/7`.

### Unit testing and continuous integration (CI)

This project is synchronized with Travis CI such that each PR gets tested before
merging (and the build is automatically triggered after each new commit).
For the maintainability of this project, it is important to understand and fix
the failing doctests if they exist.
We develop in Julia v0.6.0, but for experimentation we also build on the nightly
branch.

When you modify code in this package, you should make sure that all unit tests
pass.
To run the unit tests locally, you should do:

```
$ julia --color=yes test/runtests.jl
```

Alternatively, you can achieve the same from inside the REPL using the following
command:

```julia
julia> Pkg.test("LazySets")
```

We also advise adding new unit tests when adding new features to ensure
long-term support of your contributions.

### Contributing to the documentation

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

This documentation you are currently reading is written in Markdown, and it
relies on [Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/) to
produce the HTML layout.
The sources for creating this documentation are found in `docs/src`.
You can easily include the documentation that you wrote for your functions or
types there (see the
[Documenter.jl guide](https://juliadocs.github.io/Documenter.jl/stable/man/guide/)
or our sources for examples).

To generate the documentation locally, run `make.jl`, e.g., by executing the
following command in the terminal:

```
$ julia --color=yes docs/make.jl
```

Note that this also runs all doctests which will take some time.

## Related projects

The project [3PLIB](https://3plib.wordpress.com/) is a Java Library developed
by Frédéric Viry, and it is one of the previous works that led to the creation
of `LazySets.jl`.
3PLIB is specialized to planar projections of convex polyhedra.
It was initially created to embed this feature in Java applications, and also
provides a backend for visualization of high-dimensional reach set
approximations computed with [SpaceEx](http://spaceex.imag.fr/).

## Credits

Here we list the names of the maintainers of the `LazySets.jl` library, as well as past and present contributors (in alphabetic order).

### Core developers

- [Marcelo Forets](http://main.marcelo-forets.fr/)
- [Christian Schilling](https://schillic.github.io/)

### Contributors

- Tomer Arnon
- Kostiantyn Potomkin
- Frédéric Viry

### Colaborators

We are also grateful to the following persons for enlightening discussions: 

- [Sergiy Bogomolov](https://www.sergiybogomolov.com/)
- [Goran Frehse](https://sites.google.com/site/frehseg/) 
