```@contents
Pages = ["AbstractPolygon.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# [Polygons (AbstractPolygon)](@id def_AbstractPolygon)

A polygon is a two-dimensional polytope.

```@docs
AbstractPolygon
```

This interface requires to implement the following functions:

```@docs
tohrep(::AbstractPolygon)
tovrep(::AbstractPolygon)
```

This interface defines the following functions:

```@docs
dim(::AbstractPolygon)
volume(::AbstractPolygon)
```

The following helper functions are used for sorting directions:

```@docs
LazySets.jump2pi
⪯(::AbstractVector, ::AbstractVector)
LazySets._leq_trig(::AbstractVector{N}, ::AbstractVector{N}) where {N<:AbstractFloat}
LazySets.quadrant(::AbstractVector{N}) where {N}
```

## Implementations

* [Polygon in vertex representation (VPolygon)](@ref def_VPolygon)
