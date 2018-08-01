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

Every convex set in this library implements this interface.

```@docs
LazySet
```

### Support function and support vector

Every `LazySet` type must define a function `σ` to compute the support vector.

```@docs
support_vector
ρ(::AbstractVector{Real}, ::LazySet{Real})
support_function
```

### Other globally defined set functions

```@docs
norm(::LazySet, ::Real)
radius(::LazySet, ::Real)
diameter(::LazySet, ::Real)
an_element(::LazySet{Real})
```

## Point symmetric set

Point symmetric sets such as balls of different norms are characterized by a
center.
Note that there is a special interface combination
[Point symmetric polytope](@ref).

```@docs
AbstractPointSymmetric
```

This interface defines the following functions:

```@docs
dim(::AbstractPointSymmetric)
an_element(::AbstractPointSymmetric{N}) where {N<:Real}
```

## Polytope

A polytope has finitely many vertices (*V-representation*) resp. facets
(*H-representation*).
Note that there is a special interface combination
[Point symmetric polytope](@ref).

```@docs
AbstractPolytope
```

This interface defines the following functions:

```@docs
singleton_list(::AbstractPolytope)
linear_map(::AbstractMatrix, ::AbstractPolytope)
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

This interface defines the following functions:

```@docs
dim(::AbstractPointSymmetricPolytope)
an_element(::AbstractPointSymmetricPolytope{N}) where {N<:Real}
```

#### Hyperrectangle

A hyperrectangle is a special point symmetric polytope with axis-aligned facets.

```@docs
AbstractHyperrectangle
```

This interface defines the following functions:

```@docs
norm(::AbstractHyperrectangle, ::Real)
radius(::AbstractHyperrectangle, ::Real)
σ(::AbstractVector{Real}, ::AbstractHyperrectangle{Real})
∈(::AbstractVector{Real}, ::AbstractHyperrectangle{Real})
vertices_list(::AbstractHyperrectangle{Real})
```

#### Singleton

A singleton is a special hyperrectangle consisting of only one point.

```@docs
AbstractSingleton
```

This interface defines the following functions:

```@docs
σ(::AbstractVector{N}, ::AbstractSingleton{N}) where {N<:Real}
∈(::AbstractVector{N}, ::AbstractSingleton{N}) where {N<:Real}
an_element(::AbstractSingleton{N}) where {N<:Real}
center(::AbstractSingleton{Real})
vertices_list(::AbstractSingleton{N}) where {N<:Real}
radius_hyperrectangle(::AbstractSingleton{Real})
radius_hyperrectangle(::AbstractSingleton{Real}, ::Int)
linear_map(::AbstractMatrix, ::AbstractSingleton{N}) where {N<:Real}
```
