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
⊆(::LazySet{Float64}, ::AbstractHyperrectangle{Float64})
⊆(::AbstractPolytope{Float64}, ::AbstractHyperrectangle{Float64})
⊆(::AbstractHyperrectangle{Float64}, ::AbstractHyperrectangle{Float64})
⊆(::AbstractPolytope{Float64}, ::LazySet{Float64})
⊆(::AbstractSingleton{Float64}, ::LazySet{Float64})
⊆(::AbstractSingleton{Float64}, ::AbstractHyperrectangle{Float64})
⊆(::AbstractSingleton{Float64}, ::AbstractSingleton{Float64})
⊆(::Ball2{Float64}, ::Ball2{Float64})
⊆(::Ball2{Float64}, ::AbstractSingleton{Float64})
```

## Check for emptiness of intersection

```@docs
is_intersection_empty(::AbstractHyperrectangle{Float64}, ::AbstractHyperrectangle{Float64})
is_intersection_empty(::AbstractSingleton{Float64}, ::AbstractHyperrectangle{Float64})
is_intersection_empty(::AbstractHyperrectangle{Float64}, ::AbstractSingleton{Float64})
is_intersection_empty(::AbstractSingleton{Float64}, ::LazySet{Float64})
is_intersection_empty(::LazySet{Float64}, ::AbstractSingleton{Float64})
is_intersection_empty(::AbstractSingleton{Float64}, ::AbstractSingleton{Float64})
is_intersection_empty(::Zonotope{Float64}, ::Hyperplane{Float64})
is_intersection_empty(::Hyperplane{Float64}, ::Zonotope{Float64})
is_intersection_empty(::Ball2{Float64}, ::Ball2{Float64})
is_intersection_empty(::LineSegment{Float64}, ::LineSegment{Float64})
```

## Intersection of two sets

```@docs
intersection(::Line{Float64}, ::Line{Float64})
intersection(::Hyperrectangle{Float64}, ::Hyperrectangle{Float64})
```
