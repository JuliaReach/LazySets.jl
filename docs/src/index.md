# LazySets.jl

`LazySets` is a [Julia](http://julialang.org) package for calculus with
convex sets.

## Features

At the core of `LazySets` there are:

- Lazy (i.e. symbolic) types for several convex sets such as
  convex polygons, different classes of polytopes, and special types
  such as linear constraints.
- Most commonly used set operations, e.g. Minkowski sum, Cartesian product,
  convex hull and interval hull approximations. Moreover, lazy linear maps and
  lazy exponential maps are also provided.

Each instance of the abstract type `LazySet` implements a function,
$\sigma(d, \mathcal{X})$, to compute the supoprt vector of a set
$\mathcal{X}\subset \mathbb{R}^n$ in a given (arbitrary) direction
$d \in \mathbb{R}^n$.
This has the advantage of being able to perform only the required operations
on-demand.

On top of the previous basic type representations and operations, the following
functionality is available:

- Efficient evaluation of the support vector of nested lazy sets.
- Cartesian decomposition of lazy sets using support vectors.
- Fast overapproximation of symbolic set computations using a polyhedral
  approximation.


## Manual Outline

```@contents
Pages = [
    "man/getting_started.md",
    "man/representations.md",
    "man/operations.md",
    "man/approximations.md"
]
Depth = 2
```

