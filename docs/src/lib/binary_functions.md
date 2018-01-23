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
⊆
is_subset(::LazySet{Float64}, ::AbstractHyperrectangle{Float64})
is_subset(::AbstractPolytope{Float64}, ::AbstractHyperrectangle{Float64})
is_subset(::AbstractHyperrectangle{Float64}, ::AbstractHyperrectangle{Float64})
is_subset(::AbstractPolytope{Float64}, ::LazySet{Float64})
is_subset(::AbstractSingleton{Float64}, ::LazySet{Float64})
is_subset(::AbstractSingleton{Float64}, ::AbstractHyperrectangle{Float64})
is_subset(::AbstractSingleton{Float64}, ::AbstractSingleton{Float64})
is_subset(::Ball2{Float64}, ::Ball2{Float64})
is_subset(::Ball2{Float64}, ::AbstractSingleton{Float64})
```

## Check for emptiness of intersection

```@docs
is_intersection_empty(::AbstractHyperrectangle{Float64}, ::AbstractHyperrectangle{Float64})
is_intersection_empty(::AbstractSingleton{Float64}, ::AbstractHyperrectangle{Float64})
is_intersection_empty(::AbstractHyperrectangle{Float64}, ::AbstractSingleton{Float64})
is_intersection_empty(::AbstractSingleton{Float64}, ::LazySet{Float64})
is_intersection_empty(::LazySet{Float64}, ::AbstractSingleton{Float64})
is_intersection_empty(::AbstractSingleton{Float64}, ::AbstractSingleton{Float64})
is_intersection_empty(::Ball2{Float64}, ::Ball2{Float64})
```
