# Binary Functions on Sets

This section of the manual describes the binary functions for set types.

```@contents
Pages = ["binary_functions.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
DocTestSetup = quote
    using LazySets
end
```

## Subset check

```@docs
⊆(::LazySet{N}, ::AbstractHyperrectangle{N}, ::Bool=false) where {N<:Real}
⊆(::AbstractPolytope{N}, ::LazySet{N}, ::Bool=false) where {N<:Real}
⊆(::AbstractPolytope{N}, ::AbstractHyperrectangle, ::Bool=false) where {N<:Real}
⊆(::AbstractHyperrectangle{N}, ::AbstractHyperrectangle{N}, ::Bool=false) where {N<:Real}
⊆(::AbstractSingleton{N}, ::LazySet{N}, ::Bool=false) where {N<:Real}
⊆(::AbstractSingleton{N}, ::AbstractHyperrectangle{N}, ::Bool=false) where {N<:Real}
⊆(::AbstractSingleton{N}, ::AbstractSingleton{N}, ::Bool=false) where {N<:Real}
⊆(::Ball2{AbstractFloat}, ::Ball2{AbstractFloat}, ::Bool=false) where {N<:AbstractFloat}
	⊆(::Union{Ball2{N}, Ballp{N}}, ::AbstractSingleton{N}, ::Bool=false) where {N<:Real}
⊆(::LineSegment{N}, ::LazySet{N}, ::Bool=false) where {N<:Real}
⊆(::LineSegment{N}, ::Hyperrectangle{N}, ::Bool=false) where {N<:Real}
⊆(::Interval, ::Interval)
```

## Check for emptiness of intersection

```@docs
isdisjoint
is_intersection_empty(::AbstractHyperrectangle{N}, ::AbstractHyperrectangle{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::AbstractHyperrectangle{N}, ::AbstractSingleton{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::AbstractSingleton{N}, ::LazySet{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::AbstractSingleton{N}, ::AbstractHyperrectangle{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::AbstractSingleton{N}, ::AbstractSingleton{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::LazySet{N}, ::AbstractSingleton{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::Zonotope{N}, ::Hyperplane{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::Hyperplane{N}, ::Zonotope{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::Ball2{N}, ::Ball2{N}, ::Bool=false) where {N<:AbstractFloat}
is_intersection_empty(::LineSegment{N}, ::LineSegment{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::AbstractPolytope{N}, ::AbstractPolytope{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::AbstractSingleton{N}, ::AbstractPolytope{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::LazySet{N}, ::Union{Hyperplane{N}, Line{N}}, ::Bool=false) where {N<:Real}
is_intersection_empty(::LazySet{N}, ::HalfSpace{N}, ::Bool=false) where {N<:Real}
is_intersection_empty(::LazySet{N}, ::Union{HPolytope{N}, AbstractHPolygon{N}}) where {N<:Real}
```

## Intersection of two sets

```@docs
intersection(::Line{N}, ::Line{N}) where {N<:Real}
intersection(::Hyperrectangle{N}, ::Hyperrectangle{N}) where {N<:Real}
intersection(::AbstractHPolygon{N}, ::AbstractHPolygon{N}) where {N<:Real}
```
