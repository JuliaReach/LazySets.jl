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
cartesian_product(::HPoly{N}, ::HPoly{N}) where {N<:Real}
cartesian_product(::VPolytope{N}, ::VPolytope{N}) where N
```

## Check for emptiness of intersection

```@docs
isdisjoint
is_intersection_empty(::LazySet{N}, ::LazySet{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::AbstractHyperrectangle{N}, ::AbstractHyperrectangle{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::LazySet{N}, ::AbstractSingleton{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::AbstractHyperrectangle{N}, ::AbstractSingleton{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::AbstractSingleton{N}, ::AbstractSingleton{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::Zonotope{N}, ::Hyperplane{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::Ball2{N}, ::Ball2{N}, ::Bool=false) where {N<:AbstractFloat}
is_intersection_empty(::LineSegment{N}, ::LineSegment{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::LazySet{N}, ::Union{Hyperplane{N}, Line2D{N}}, ::Bool=false) where {N<:Real}
is_intersection_empty(::LazySet{N}, ::HalfSpace{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::HalfSpace{N}, ::HalfSpace{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::AbstractPolyhedron{N}, ::LazySet{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::UnionSet{N}, ::LazySet{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::UnionSetArray{N}, ::LazySet{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::Universe{N}, ::LazySet{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::Complement{N}, ::LazySet{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::Zonotope{N}, ::Zonotope{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::Interval{N}, ::Interval{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::CartesianProductArray{N}, ::AbstractPolyhedron{N}) where {N<:Real}
is_intersection_empty(::CartesianProductArray{N}, ::CartesianProductArray{N}) where {N<:Real}
is_intersection_empty(::CartesianProductArray{N}, ::AbstractHyperrectangle{N}, ::Bool=false) where {N<:Real}
```

## Convex hull

```@docs
convex_hull(::LazySet{N}, ::LazySet{N}) where {N<:Real}
convex_hull(::HPoly{N}, ::HPoly{N}) where {N<:Real}
convex_hull(::VPolytope{N}, ::VPolytope{N}) where {N<:Real}
convex_hull(::VPolygon{N}, ::VPolygon{N}) where {N<:Real}
convex_hull(::Vector{VN}) where {N<:Real, VN<:AbstractVector{N}}
convex_hull(U::UnionSetArray{N, PT}) where {N, PT<:AbstractPolytope{N}}
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
intersection(::Union{VPolytope{N}, VPolygon{N}}, ::Union{VPolytope{N}, VPolygon{N}}) where {N<:Real}
intersection(::UnionSet{N}, ::LazySet{N}) where {N<:Real}
intersection(::UnionSetArray{N}, ::LazySet{N}) where {N<:Real}
intersection(::Universe{N}, ::LazySet{N}) where {N<:Real}
intersection(::AbstractPolyhedron{N}, ::ResetMap{N}) where {N<:Real}
intersection(::CartesianProductArray{N}, ::CartesianProductArray{N}) where {N<:Real}
intersection(::LinearMap{N}, ::LazySet{N}) where {N<:Real}
intersection(::CartesianProductArray{N}, ::AbstractPolyhedron{N}) where {N<:Real}
```

## Minkowski sum

```@docs
minkowski_sum(::LazySet{N}, ::LazySet{N}) where {N<:Real}
minkowski_sum(::AbstractPolyhedron{N}, ::AbstractPolyhedron{N}) where {N<:Real}
minkowski_sum(::VPolytope{N}, ::VPolytope{N}) where {N<:Real}
minkowski_sum(::AbstractHyperrectangle{N}, ::AbstractHyperrectangle{N}) where {N<:Real}
minkowski_sum(::AbstractZonotope{N}, ::AbstractZonotope{N}) where {N<:Real}
minkowski_sum(::VPolygon{N}, ::VPolygon{N}) where {N<:Real}
minkowski_sum(::PolynomialZonotope, ::Zonotope)
minkowski_sum(::Interval, ::Interval)
```

## Minkowski difference
```@docs
minkowski_difference(::LazySet{N}, ::LazySet{N}) where {N<:Real}
pontryagin_difference
```

## Subset check

```@docs
⊆(::LazySet{N}, ::LazySet{N}, ::Bool=false) where {N<:Real}
⊆(::LazySet{N}, ::AbstractHyperrectangle{N}, ::Bool=false) where {N<:Real}
⊆(::AbstractPolytope{N}, ::LazySet{N}, ::Bool=false) where {N<:Real}
⊆(::AbstractPolytope{N}, ::AbstractHyperrectangle, ::Bool=false) where {N<:Real}
⊆(::AbstractHyperrectangle{N}, ::AbstractHyperrectangle{N}, ::Bool=false) where {N<:Real}
⊆(::LazySet{N}, ::AbstractPolyhedron{N}, ::Bool=false) where {N<:Real}
⊆(::AbstractSingleton{N}, ::LazySet{N}, ::Bool=false) where {N<:Real}
⊆(::AbstractSingleton{N}, ::AbstractHyperrectangle{N}, ::Bool=false) where {N<:Real}
⊆(::AbstractSingleton{N}, ::AbstractSingleton{N}, ::Bool=false) where {N<:Real}
⊆(::Ball2{N}, ::Ball2{N}, ::Bool=false) where {N<:AbstractFloat}
⊆(::Union{Ball2{N}, Ballp{N}}, ::AbstractSingleton{N}, ::Bool=false) where {N<:AbstractFloat}
⊆(::LineSegment{N}, ::LazySet{N}, ::Bool=false) where {N<:Real}
⊆(::LineSegment{N}, ::AbstractHyperrectangle{N}, ::Bool=false) where {N<:Real}
⊆(::Interval, ::Interval)
⊆(::EmptySet{N}, ::LazySet{N}, ::Bool=false) where {N<:Real}
⊆(::LazySet{N}, ::EmptySet{N}, ::Bool=false) where {N<:Real}
⊆(::UnionSet{N}, ::LazySet{N}, ::Bool=false) where {N<:Real}
⊆(::UnionSetArray{N}, ::LazySet{N}, ::Bool=false) where {N<:Real}
⊆(::LazySet{N}, ::Universe{N}, ::Bool=false) where {N<:Real}
⊆(::Universe{N}, ::LazySet{N}, ::Bool=false) where {N<:Real}
⊆(::LazySet{N}, ::Complement{N}, ::Bool=false) where {N<:Real}
⊆(::CartesianProductArray{N}, ::CartesianProductArray{N}, ::Bool=false) where {N<:Real}
```

## Set difference

```@docs
\(::LazySet, ::LazySet)
difference(::IN, ::IN) where {N, IN<:Interval{N}}
difference(::AbstractHyperrectangle{N}, ::AbstractHyperrectangle{N}) where {N}
```
