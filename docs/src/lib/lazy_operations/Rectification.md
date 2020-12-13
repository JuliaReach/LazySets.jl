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
σ(::AbstractVector, ::Rectification)
σ(::AbstractVector, ::Rectification{N, <:AbstractHyperrectangle{N}}) where {N}
σ(::AbstractVector, ::Rectification{N, <:CartesianProduct{N}}) where {N}
σ(::AbstractVector, ::Rectification{N, <:CartesianProductArray{N}}) where {N}
ρ(::AbstractVector, ::Rectification)
an_element(::Rectification)
∈(::AbstractVector, ::Rectification)
isempty(::Rectification)
isbounded(::Rectification{N}) where {N}
to_union_of_projections(::Rectification, ::Bool=false)
```

Inherited from [`LazySet`](@ref):
* [`singleton_list`](@ref singleton_list(::LazySet))

## Rectification cache

```@docs
RectificationCache
```
