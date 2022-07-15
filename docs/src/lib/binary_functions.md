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
cartesian_product(::ConvexSet, ::ConvexSet)
```

## Check for emptiness of intersection

!!! note
    `is_intersection_empty` can be used as an alternative name to `isdisjoint`.

```@docs
isdisjoint(::ConvexSet, ::ConvexSet, ::Bool=false)
isdisjoint(::AbstractHyperrectangle, ::AbstractHyperrectangle, ::Bool=false)
isdisjoint(::ConvexSet, ::AbstractSingleton, ::Bool=false)
isdisjoint(::AbstractHyperrectangle, ::AbstractSingleton, ::Bool=false)
isdisjoint(::AbstractSingleton, ::AbstractSingleton, ::Bool=false)
isdisjoint(::AbstractZonotope, ::Union{Hyperplane, Line2D}, ::Bool=false)
isdisjoint(::Ball2, ::Ball2, ::Bool=false)
isdisjoint(::LineSegment, ::LineSegment, ::Bool=false)
isdisjoint(::ConvexSet, ::Union{Hyperplane, Line2D}, ::Bool=false)
isdisjoint(::ConvexSet, ::HalfSpace, ::Bool=false)
isdisjoint(::HalfSpace, ::HalfSpace, ::Bool=false)
isdisjoint(::AbstractPolyhedron, ::ConvexSet, ::Bool=false)
isdisjoint(::UnionSet, ::ConvexSet, ::Bool=false)
isdisjoint(::UnionSetArray, ::ConvexSet, ::Bool=false)
isdisjoint(::Universe, ::ConvexSet, ::Bool=false)
isdisjoint(::Complement, ::ConvexSet, ::Bool=false)
isdisjoint(::AbstractZonotope, ::AbstractZonotope, ::Bool=false)
isdisjoint(::Interval, ::Interval, ::Bool=false)
isdisjoint(::CartesianProductArray, ::AbstractPolyhedron)
isdisjoint(::CartesianProductArray, ::CartesianProductArray)
isdisjoint(::CartesianProductArray, ::AbstractHyperrectangle, ::Bool=false)
isdisjoint(::Line2D, ::Line2D, ::Bool=false)
```

## Convex hull

```@docs
convex_hull(::ConvexSet{N}, ::ConvexSet{N}) where {N}
convex_hull(::HPoly, ::HPoly)
convex_hull(::VPolytope, ::VPolytope)
convex_hull(::VPolygon, ::VPolygon)
convex_hull(::Vector{VN}) where {N, VN<:AbstractVector{N}}
convex_hull(::UnionSetArray{N, PT}; kwargs...) where {N, PT<:AbstractPolytope{N}}
monotone_chain!
```

## Intersection of two sets

```@docs
intersection(::AbstractSingleton, ::ConvexSet)
intersection(::Line2D, ::Line2D)
intersection(::AbstractHyperrectangle, ::AbstractHyperrectangle)
intersection(::Interval, ::Interval)
intersection(::Interval, ::HalfSpace)
intersection(::Interval, ::Hyperplane)
intersection(::Interval, ::ConvexSet)
intersection(::AbstractHPolygon, ::AbstractHPolygon)
intersection(::AbstractPolyhedron{N}, ::AbstractPolyhedron{N}) where {N}
intersection(::Union{VPolytope, VPolygon}, ::Union{VPolytope, VPolygon})
intersection(::VPolygon, ::VPolygon; ::Bool=true)
intersection(::UnionSet, ::ConvexSet)
intersection(::UnionSetArray, ::ConvexSet)
intersection(::Universe, ::ConvexSet)
intersection(::AbstractPolyhedron, ::ResetMap)
intersection(::CartesianProductArray, ::CartesianProductArray)
intersection(::LinearMap, ::ConvexSet)
intersection(::CartesianProductArray, ::AbstractPolyhedron)
intersection(::LineSegment, ::Line2D)
intersection(::LineSegment, ::LineSegment)
intersection(::AbstractZonotope{N}, ::HalfSpace{N}) where {N}
```

## Minkowski sum

```@docs
minkowski_sum(::ConvexSet, ::ConvexSet)
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
minkowski_difference(::ConvexSet, ::ConvexSet)
pontryagin_difference
```

## Subset check

```@docs
issubset
⊆(::ConvexSet, ::ConvexSet, ::Bool=false)
⊆(::ConvexSet, ::AbstractHyperrectangle, ::Bool=false)
⊆(::AbstractPolytope, ::ConvexSet, ::Bool=false)
⊆(::AbstractPolytope, ::AbstractHyperrectangle, ::Bool=false)
⊆(::AbstractZonotope, ::AbstractHyperrectangle)
⊆(::AbstractHyperrectangle, ::AbstractHyperrectangle, ::Bool=false)
⊆(::ConvexSet, ::AbstractPolyhedron, ::Bool=false)
⊆(::AbstractSingleton, ::ConvexSet, ::Bool=false)
⊆(::AbstractSingleton, ::AbstractHyperrectangle, ::Bool=false)
⊆(::AbstractSingleton, ::AbstractSingleton, ::Bool=false)
⊆(::Ball2, ::Ball2, ::Bool=false)
⊆(::Union{Ball2, Ballp}, ::AbstractSingleton, ::Bool=false)
⊆(::LineSegment, ::ConvexSet, ::Bool=false)
⊆(::LineSegment, ::AbstractHyperrectangle, ::Bool=false)
⊆(::Interval{N}, ::Interval, ::Bool=false) where {N}
⊆(::Interval, ::UnionSet, ::Bool=false)
⊆(::EmptySet, ::ConvexSet, ::Bool=false)
⊆(::ConvexSet, ::EmptySet, ::Bool=false)
⊆(::UnionSet, ::ConvexSet, ::Bool=false)
⊆(::UnionSetArray, ::ConvexSet, ::Bool=false)
⊆(::ConvexSet, ::Universe, ::Bool=false)
⊆(::Universe, ::ConvexSet, ::Bool=false)
⊆(::ConvexSet, ::Complement, ::Bool=false)
⊆(::CartesianProduct, ::CartesianProduct, ::Bool=false)
⊆(::CartesianProductArray, ::CartesianProductArray, ::Bool=false)
⊆(::AbstractZonotope, ::AbstractHyperrectangle, ::Bool=false)
⊂
```

## Set difference

```@docs
\(::ConvexSet, ::ConvexSet)
difference(::IN, ::IN) where {N, IN<:Interval{N}}
difference(::AbstractHyperrectangle{N}, ::AbstractHyperrectangle{N}) where {N}
```

## Distance

```@docs
distance(::AbstractSingleton, ::ConvexSet; ::Real=2.0)
distance(::AbstractHyperrectangle, ::AbstractHyperrectangle; ::Real=2.0)
```
