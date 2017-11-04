# LazySets.jl

`LazySets` is a [Julia](http://julialang.org) package for calculus with
convex sets.

The aim is to provide a scalable library for solving complex set-based
problems, such as those encountered in [differential inclusions](https://en.wikipedia.org/wiki/Differential_inclusion)
or reachability analysis techniques in the domain of [formal verification](https://en.wikipedia.org/wiki/Formal_verification).
Typically, one has to solve a set-based recurrence, and for visualization
purposes the final result is obtained through an adequate projection onto low-dimensions.
This library implements types to construct the set formulas with the usual mathematical
notation, and methods to efficiently and accurately approximate the
projection in low-dimensions. 

The strategy consists of using lazy (i.e. symbolic) representations of complex
set formulas. This provides an exact but abstract representation of convex sets and common
operations. Then, concrete information is obtained through querying specific directions.
More precisely, each concrete subtype $\mathcal{X}$ of the abstract type `LazySet`,
exports a method to calculate its support vector $\sigma(d, \mathcal{X})$ in a given (arbitrary) direction
$d \in \mathbb{R}^n$. Representing sets exactly but lazily has the advantage of being
able to perform only the required operations on-demand.

For very long computations (e.g. set-based recurrences with tens of thousands of elements),
it is useful to combine both lazy and concrete representations such as polyhedral approximations.
All this is easy to do with `LazySets`. Moreover, there is a specialized module for handling
two-dimensional projections using Cartesian decomposition techniques. The projection
can be taken to the desired precision using an iterative refinement method.

## Example



## Features

At the core of `LazySets` there is the following functionality:

- Lazy (i.e. symbolic) types for several classes of convex sets such as
  balls in different norms, polygons in constraint or vertex representation,
  special types such as lines and linear constraints, hyperrectangles, and
  high-dimensional polyhedra.
- Most commonly used set operations, e.g. Minkowski sum, Cartesian product,
  convex hull and interval hull approximations. Moreover, lazy linear maps and
  lazy exponential maps are also provided.

On top of the previous basic type representations and operations, `LazySets`
can be used to:

- Efficiently evaluate the support vector of nested lazy sets using parametrized
`LazySet` arrays.
- Cartesian decomposition of lazy sets using two-dimensional projections.
- Fast overapproximation of an exact set using a polyhedral
  approximation, to the desired accuracy.

## Manual Outline

```@contents
Pages = [
    "man/getting_started.md",
    "man/support_vectors.md",
    "man/polyhedral_approximations.md",
    "man/fast_2d_LPs.md",
    "man/iterative_refinement.md"
]
Depth = 2
```

## Library Outline

```@contents
Pages = [
    "lib/representations.md",
    "lib/operations.md",
    "lib/approximations.md",
    "lib/utils.md"
]
Depth = 2
```
