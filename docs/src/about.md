# About

This page contains some general information about this project, and
recommendations about contributing.

```@contents
Pages = ["about.md"]
```

## Contributing

If you like this package, consider contributing!

[Creating an issue](https://help.github.com/en/articles/creating-an-issue) in the [LazySets GitHub issue tracker](https://github.com/JuliaReach/LazySets.jl/issues) to report a bug, open a discussion about existing functionality or suggesting new one is appreciated.

If you have written code and would like it to be peer reviewed and added to the library, you can [fork](https://help.github.com/en/articles/fork-a-repo) the repository and send a pull request (see below). Typical contributions include fixing a bug, adding a new feature or improving the documentation (either in source code or the online manual).

You are also welcome to get in touch with us in the [JuliaReach gitter chat](https://gitter.im/JuliaReach/Lobby).

Below we detail some general comments about contributing to this package. The [JuliaReach Developer's Documentation](https://juliareach.github.io/JuliaReachDevDocs/latest/) describes coding guidelines; take a look when in doubt about the coding style that is expected for the code that is finally merged into the library.

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

New functions and types should be documented according to the
[JuliaReach Developer's Documentation](https://juliareach.github.io/JuliaReachDevDocs/latest/guidelines/#Writing-docstrings-1).

You can view the source code documentation from inside the REPL by typing `?`
followed by the name of the type or function.
For example, the following command will print the documentation of the `LazySet`
type:

```julia
julia> ?LazySet
```

The documentation you are currently reading is written in [Markdown](https://en.wikipedia.org/wiki/Markdown), and it
relies on the package [Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/) to
produce the HTML layout.
The sources for creating this documentation are found in `docs/src`.
You can easily include the documentation that you wrote for your functions or
types there (see the source code or [Documenter's guide](https://juliadocs.github.io/Documenter.jl/stable/man/guide/)for examples).

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
- Sebastián Guadalupe
- Kostiantyn Potomkin
- Frédéric Viry

### Colaborators

We are also grateful to the following persons for enlightening discussions: 

- [Sergiy Bogomolov](https://www.sergiybogomolov.com/)
- [Goran Frehse](https://sites.google.com/site/frehseg/) 
