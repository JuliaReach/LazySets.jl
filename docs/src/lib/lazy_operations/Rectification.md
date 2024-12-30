```@meta
CurrentModule = LazySets
```

# [Rectification](@id def_Rectification)

```@docs
Rectification
set(::Rectification)
dim(::Rectification)
σ(::AbstractVector, ::Rectification)
σ(::AbstractVector, ::Rectification{N, <:AbstractHyperrectangle}) where {N}
σ(::AbstractVector, ::Rectification{N, <:CartesianProduct}) where {N}
σ(::AbstractVector, ::Rectification{N, <:CartesianProductArray}) where {N}
ρ(::AbstractVector, ::Rectification)
an_element(::Rectification)
∈(::AbstractVector, ::Rectification)
isempty(::Rectification)
isbounded(::Rectification)
to_union_of_projections(::Rectification, ::Bool=false)
```
Inherited from [`LazySet`](@ref):
* [`singleton_list`](@ref singleton_list(::LazySet))

## Rectification cache

```@docs
RectificationCache
```
