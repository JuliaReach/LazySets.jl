# Binary Functions on Sets

This section of the manual describes the binary functions for set types.

```@contents
Pages = ["binary_functions.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

## Cartesian product

```@docs
cartesian_product(::VPolytope, ::VPolytope)
cartesian_product(::LazySet, ::LazySet)
```

## Check for emptiness of intersection

!!! note
    `is_intersection_empty` can be used as an alternative name to `isdisjoint`.

```@docs
isdisjoint(::LazySet, ::LazySet, ::Bool=false)
isdisjoint(::AbstractHyperrectangle, ::AbstractHyperrectangle, ::Bool=false)
isdisjoint(::LazySet, ::AbstractSingleton, ::Bool=false)
isdisjoint(::AbstractHyperrectangle, ::AbstractSingleton, ::Bool=false)
isdisjoint(::AbstractSingleton, ::AbstractSingleton, ::Bool=false)
isdisjoint(::AbstractZonotope, ::Union{Hyperplane, Line2D}, ::Bool=false)
isdisjoint(::Ball2, ::Ball2, ::Bool=false)
isdisjoint(::LineSegment, ::LineSegment, ::Bool=false)
isdisjoint(::LazySet, ::Union{Hyperplane, Line2D}, ::Bool=false)
isdisjoint(::LazySet, ::HalfSpace, ::Bool=false)
isdisjoint(::HalfSpace, ::HalfSpace, ::Bool=false)
isdisjoint(::AbstractPolyhedron, ::LazySet, ::Bool=false)
isdisjoint(::UnionSet, ::LazySet, ::Bool=false)
isdisjoint(::UnionSetArray, ::LazySet, ::Bool=false)
isdisjoint(::Universe, ::LazySet, ::Bool=false)
isdisjoint(::Complement, ::LazySet, ::Bool=false)
isdisjoint(::AbstractZonotope, ::AbstractZonotope, ::Bool=false)
isdisjoint(::Interval, ::Interval, ::Bool=false)
isdisjoint(::CartesianProductArray, ::AbstractPolyhedron)
isdisjoint(::CartesianProductArray, ::CartesianProductArray)
isdisjoint(::CartesianProductArray, ::AbstractHyperrectangle, ::Bool=false)
isdisjoint(::Line2D, ::Line2D, ::Bool=false)
```

## Convex hull

```@docs
convex_hull(::LazySet{N}, ::LazySet{N}) where {N}
convex_hull(::HPoly, ::HPoly)
convex_hull(::VPolytope, ::VPolytope)
convex_hull(::VPolygon, ::VPolygon)
convex_hull(::Vector{VN}) where {N, VN<:AbstractVector{N}}
convex_hull(::UnionSetArray{N, PT}; kwargs...) where {N, PT<:AbstractPolytope{N}}
monotone_chain!
```

## Intersection of two sets

```@docs
intersection(::AbstractSingleton, ::LazySet)
intersection(::Line2D, ::Line2D)
intersection(::AbstractHyperrectangle, ::AbstractHyperrectangle)
intersection(::Interval, ::Interval)
intersection(::Interval, ::HalfSpace)
intersection(::Interval, ::Hyperplane)
intersection(::Interval, ::LazySet)
intersection(::AbstractHPolygon, ::AbstractHPolygon)
intersection(::AbstractPolyhedron{N}, ::AbstractPolyhedron{N}) where {N}
intersection(::Union{VPolytope, VPolygon}, ::Union{VPolytope, VPolygon})
intersection(::VPolygon, ::VPolygon; ::Bool=true)
intersection(::UnionSet, ::LazySet)
intersection(::UnionSetArray, ::LazySet)
intersection(::Universe, ::LazySet)
intersection(::AbstractPolyhedron, ::ResetMap)
intersection(::CartesianProductArray, ::CartesianProductArray)
intersection(::LinearMap, ::LazySet)
intersection(::CartesianProductArray, ::AbstractPolyhedron)
intersection(::LineSegment, ::Line2D)
intersection(::LineSegment, ::LineSegment)
intersection(::AbstractZonotope{N}, ::HalfSpace{N}) where {N}
```

## Minkowski sum

```@docs
minkowski_sum(::LazySet, ::LazySet)
minkowski_sum(::AbstractPolyhedron, ::AbstractPolyhedron)
minkowski_sum(::VPolytope, ::VPolytope)
minkowski_sum(::AbstractHyperrectangle, ::AbstractHyperrectangle)
minkowski_sum(::AbstractZonotope, ::AbstractZonotope)
minkowski_sum(::VPolygon, ::VPolygon)
minkowski_sum(::DensePolynomialZonotope, ::AbstractZonotope)
minkowski_sum(::Interval, ::Interval)
minkowski_sum(::AbstractSingleton, ::AbstractSingleton)
```

## Minkowski difference
```@docs
minkowski_difference(::LazySet, ::LazySet)
pontryagin_difference
```

## Subset check

```@docs
issubset
⊆(::LazySet, ::LazySet, ::Bool=false)
⊆(::LazySet, ::AbstractHyperrectangle, ::Bool=false)
⊆(::AbstractPolytope, ::LazySet, ::Bool=false)
⊆(::AbstractPolytope, ::AbstractHyperrectangle, ::Bool=false)
⊆(::AbstractZonotope, ::AbstractHyperrectangle)
⊆(::AbstractHyperrectangle, ::AbstractHyperrectangle, ::Bool=false)
⊆(::LazySet, ::AbstractPolyhedron, ::Bool=false)
⊆(::AbstractSingleton, ::LazySet, ::Bool=false)
⊆(::AbstractSingleton, ::AbstractHyperrectangle, ::Bool=false)
⊆(::AbstractSingleton, ::AbstractSingleton, ::Bool=false)
⊆(::Ball2, ::Ball2, ::Bool=false)
⊆(::Union{Ball2, Ballp}, ::AbstractSingleton, ::Bool=false)
⊆(::LineSegment, ::LazySet, ::Bool=false)
⊆(::LineSegment, ::AbstractHyperrectangle, ::Bool=false)
⊆(::Interval{N}, ::Interval, ::Bool=false) where {N}
⊆(::Interval, ::UnionSet, ::Bool=false)
⊆(::EmptySet, ::LazySet, ::Bool=false)
⊆(::LazySet, ::EmptySet, ::Bool=false)
⊆(::UnionSet, ::LazySet, ::Bool=false)
⊆(::UnionSetArray, ::LazySet, ::Bool=false)
⊆(::LazySet, ::Universe, ::Bool=false)
⊆(::Universe, ::LazySet, ::Bool=false)
⊆(::LazySet, ::Complement, ::Bool=false)
⊆(::CartesianProduct, ::CartesianProduct, ::Bool=false)
⊆(::CartesianProductArray, ::CartesianProductArray, ::Bool=false)
⊆(::AbstractZonotope, ::AbstractHyperrectangle, ::Bool=false)
⊂
```

## Set difference

```@docs
\(::LazySet, ::LazySet)
difference(::IN, ::IN) where {N, IN<:Interval{N}}
difference(::AbstractHyperrectangle{N}, ::AbstractHyperrectangle{N}) where {N}
```

## Distance

```@docs
distance(::AbstractSingleton, ::LazySet; ::Real=2.0)
distance(::AbstractHyperrectangle, ::AbstractHyperrectangle; ::Real=2.0)
```
