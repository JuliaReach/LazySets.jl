# LazySets.jl

`LazySets` is a [Julia](http://julialang.org) package for calculus with convex sets.


## Features

The core of this package are the following types:

- Lazy (i.e. symbolic) types for most commonly used convex sets, such as:
    - Convex Polygons
    - Polyhedra, e.g. Unit balls in different norms, Hyperrectangles
- Most commonly used set operations, such as:
    - Minkowski sum
    - Cartesian product
    - Convex hull
    - Linear maps and exponential maps

On top of this, the following functionality is provided:

- Efficient evaluation of the support vector of nested lazy sets 
- Cartesian decomposition of a convex set
- Efficient projection from high-dimensions to two-dimensions using Lotov's
  algorithm

In `LazySets`, each set and operation provides a function to compute its support
vector in a given (arbitrary) direction. This has the advantage of being able to 
perform only the required operations on-demand.

In the rest of this section we recall the required preliminaries from convex analysis. 
 
## Support function and support vector



