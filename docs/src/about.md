# Contributing

```@contents
Pages = ["about/CONTRIBUTING.md"]
```

This page details the some of the guidelines that should be followed when contributing to this package.

## Running the Unit Tests

```julia
$ julia --color=yes test/runtests.jl
```

## Branches



## Contributing to the Documentation

The documentation source is written with Markdown, and we use
[Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) to produce the HTML
documentation. To build the docs, run `make.jl`:

```julia
$ julia --color=yes docs/make.jl
```

## Credits

These persons have contributed to `LazySets.jl` (in alphabetic order):

- [Marcelo Forets](http://marcelo-forets.fr)
- Christian Schilling
- Frederic Viry

We are also grateful to Goran Frehse for enlightening discussions.

## Linked Projects

The following projects are linked to `LazySets.jl`:
- [3PLIB](https://3plib.wordpress.com/) is a Java Library developed by Frederic Viry, specialized in the planar projections of convex polyhedra. It has initially been created to embed this need in Java applications, such as the results viewer of [SpaceEx](http://spaceex.imag.fr/). It is one of the previous works that led to the creation of `LazySets.jl`.
