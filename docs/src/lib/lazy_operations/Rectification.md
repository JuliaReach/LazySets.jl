```@meta
CurrentModule = LazySets
```

# [Rectification](@id def_Rectification)

Note that the rectification of a convex set is generally not convex.
Hence this set type is not part of the convex-set family `LazySet`.

```@docs
Rectification
set(::Rectification)
dim(::Rectification)
σ(::AbstractVector{N}, ::Rectification{N}) where {N<:Real}
σ(::AbstractVector{N}, ::Rectification{N, <:AbstractHyperrectangle{N}}) where {N<:Real}
σ(::AbstractVector{N}, ::Rectification{N, <:CartesianProduct{N}}) where {N<:Real}
σ(::AbstractVector{N}, ::Rectification{N, <:CartesianProductArray{N}}) where {N<:Real}
ρ(::AbstractVector{N}, ::Rectification{N}) where {N<:Real}
an_element(::Rectification{N}) where {N<:Real}
∈(::AbstractVector{N}, ::Rectification{N}) where {N<:Real}
isempty(::Rectification)
isbounded(::Rectification{N}) where {N<:Real}
to_union_of_projections(::Rectification{N}, ::Bool=false) where {N<:Real}
```

## Rectification cache

```@docs
RectificationCache
```
