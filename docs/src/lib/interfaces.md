# Set Interfaces

This section of the manual describes the interfaces for different set types.
Every set that fits the description of an interface should also implement it.
This helps in several ways:
- avoid code duplicates,
- provide functions for many sets at once,
- allow changes in the source code without changing the API.

The interface functions are outlined in the interface documentation.
See [Common Set Representations](@ref) for implementations of the interfaces.

!!! note
    The naming convention is such that all interface names (with the exception
    of the main abstract type `LazySet`) should be preceded by `Abstract`.

The following diagram shows the interface hierarchy.

![../assets/interfaces.png](../assets/interfaces.png)

```@contents
Pages = ["interfaces.md"]
Depth = 4
```

```@meta
CurrentModule = LazySets
DocTestSetup = quote
    using LazySets
end
```

## LazySet

Every set in this library implements this interface.

```@docs
LazySet
```

### Support function and support vector

Every `LazySet` type must define a function `σ` to compute the support vector.

```@docs
support_vector
ρ
support_function
```

### Other globally defined set functions

```@docs
isconvex(S::LazySet{Float64})
an_element(S::LazySet{Float64})
```

## Convex set

Common interface for basic convex sets.

```@docs
AbstractConvexSet
isconvex(S::AbstractConvexSet{Float64})
```

## Point symmetric set

Point symmetric sets such as balls of different norms are characterized by a
center.
Note that there is a special interface combination
[Point symmetric polytope](@ref).

```@docs
AbstractPointSymmetric
```

## Polytope

A polytope has finitely many vertices (*V-representation*) resp. facets
(*H-representation*).
Note that there is a special interface combination
[Point symmetric polytope](@ref).

```@docs
AbstractPolytope
```

### Polygon

A polygon is a two-dimensional polytope.

```@docs
AbstractPolygon
```

#### HPolygon

An HPolygon is a polygon in H-representation (or constraint representation).

```@docs
AbstractHPolygon
```

### Point symmetric polytope

A point symmetric polytope is a combination of two other interfaces:
[Point symmetric set](@ref) and [Polytope](@ref).

```@docs
AbstractPointSymmetricPolytope
```

#### Hyperrectangle

A hyperrectangle is a special point symmetric polytope with axis-aligned facets.

```@docs
AbstractHyperrectangle
```

#### Singleton

A singleton is a special hyperrectangle consisting of only one point.

```@docs
AbstractSingleton
```
