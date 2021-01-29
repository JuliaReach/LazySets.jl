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
cartesian_product(::HPoly, ::HPoly)
cartesian_product(::VPolytope, ::VPolytope)
cartesian_product(::LazySet, ::LazySet)
```

## Check for emptiness of intersection

!!! note
    `isdisjoint` can be used as an alternative name to `is_intersection_empty`.

```@docs
is_intersection_empty(::LazySet, ::LazySet, ::Bool=false)
is_intersection_empty(::AbstractHyperrectangle, ::AbstractHyperrectangle, ::Bool=false)
is_intersection_empty(::LazySet, ::AbstractSingleton, ::Bool=false)
is_intersection_empty(::AbstractHyperrectangle, ::AbstractSingleton, ::Bool=false)
is_intersection_empty(::AbstractSingleton, ::AbstractSingleton, ::Bool=false)
is_intersection_empty(::AbstractZonotope, ::Union{Hyperplane, Line2D}, ::Bool=false)
is_intersection_empty(::Ball2, ::Ball2, ::Bool=false)
is_intersection_empty(::LineSegment, ::LineSegment, ::Bool=false)
is_intersection_empty(::LazySet, ::Union{Hyperplane, Line2D}, ::Bool=false)
is_intersection_empty(::LazySet, ::HalfSpace, ::Bool=false)
is_intersection_empty(::HalfSpace, ::HalfSpace, ::Bool=false)
is_intersection_empty(::AbstractPolyhedron, ::LazySet, ::Bool=false)
is_intersection_empty(::UnionSet, ::LazySet, ::Bool=false)
is_intersection_empty(::UnionSetArray, ::LazySet, ::Bool=false)
is_intersection_empty(::Universe, ::LazySet, ::Bool=false)
is_intersection_empty(::Complement, ::LazySet, ::Bool=false)
is_intersection_empty(::AbstractZonotope, ::AbstractZonotope, ::Bool=false)
is_intersection_empty(::Interval, ::Interval, ::Bool=false)
is_intersection_empty(::CartesianProductArray, ::AbstractPolyhedron)
is_intersection_empty(::CartesianProductArray, ::CartesianProductArray)
is_intersection_empty(::CartesianProductArray, ::AbstractHyperrectangle, ::Bool=false)
is_intersection_empty(::Line2D, ::Line2D, ::Bool=false)
```

## Convex hull

```@docs
convex_hull(::LazySet{N}, ::LazySet{N}) where {N<:Real}
convex_hull(::HPoly, ::HPoly)
convex_hull(::VPolytope, ::VPolytope)
convex_hull(::VPolygon, ::VPolygon)
convex_hull(::Vector{VN}) where {N<:Real, VN<:AbstractVector{N}}
convex_hull(::UnionSetArray{N, PT}; kwargs...) where {N, PT<:AbstractPolytope{N}}
monotone_chain!
```

## Intersection of two sets

```@docs
intersection(::AbstractSingleton{N}, ::LazySet{N}) where {N<:Real}
intersection(::Line2D{N}, ::Line2D{N}) where {N<:Real}
intersection(::AbstractHyperrectangle{N}, ::AbstractHyperrectangle{N}) where {N<:Real}
intersection(::Interval{N}, ::Interval{N}) where {N<:Real}
intersection(::Interval{N}, ::HalfSpace{N}) where {N<:Real}
intersection(::Interval{N}, ::Hyperplane{N}) where {N<:Real}
intersection(::Interval{N}, ::LazySet{N}) where {N<:Real}
intersection(::AbstractHPolygon{N}, ::AbstractHPolygon{N}, ::Bool=true) where {N<:Real}
intersection(::AbstractPolyhedron{N}, ::AbstractPolyhedron{N}) where {N<:Real}
intersection(::Union{VPolytope, VPolygon}, ::Union{VPolytope, VPolygon})
intersection(::VPolygon{N}, ::VPolygon{N}; ::Bool=true) where {N}
intersection(::UnionSet{N}, ::LazySet{N}) where {N<:Real}
intersection(::UnionSetArray{N}, ::LazySet{N}) where {N<:Real}
intersection(::Universe{N}, ::LazySet{N}) where {N<:Real}
intersection(::AbstractPolyhedron{N}, ::ResetMap{N}) where {N<:Real}
intersection(::CartesianProductArray{N}, ::CartesianProductArray{N}) where {N<:Real}
intersection(::LinearMap{N}, ::LazySet{N}) where {N<:Real}
intersection(::CartesianProductArray{N}, ::AbstractPolyhedron{N}) where {N<:Real}
intersection(::LineSegment{N}, ::Line2D{N}) where {N<:Real}
intersection(::LineSegment{N}, ::LineSegment{N}) where {N<:Real}
```

## Minkowski sum

```@docs
minkowski_sum(::LazySet, ::LazySet)
minkowski_sum(::AbstractPolyhedron, ::AbstractPolyhedron)
minkowski_sum(::VPolytope, ::VPolytope)
minkowski_sum(::AbstractHyperrectangle, ::AbstractHyperrectangle)
minkowski_sum(::AbstractZonotope, ::AbstractZonotope)
minkowski_sum(::VPolygon, ::VPolygon)
minkowski_sum(::PolynomialZonotope, ::AbstractZonotope)
minkowski_sum(::Interval, ::Interval)
minkowski_sum(::AbstractSingleton, ::AbstractSingleton)
```

## Minkowski difference
```@docs
minkowski_difference(::LazySet{N}, ::LazySet{N}) where {N<:Real}
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
⊆(::Interval, ::Interval, ::Bool=false)
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
```

## Set difference

```@docs
\(::LazySet, ::LazySet)
difference(::IN, ::IN) where {N, IN<:Interval{N}}
difference(::AbstractHyperrectangle{N}, ::AbstractHyperrectangle{N}) where {N}
```
