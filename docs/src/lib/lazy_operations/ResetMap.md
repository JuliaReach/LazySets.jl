```@meta
CurrentModule = LazySets
```

# [Reset map (ResetMap)](@id def_ResetMap)

```@docs
ResetMap
dim(::ResetMap)
ρ(::AbstractVector{N}, ::ResetMap{N}) where {N<:Real}
σ(::AbstractVector{N}, ::ResetMap{N}) where {N<:Real}
an_element(::ResetMap)
isempty(::ResetMap)
get_A(::ResetMap{N}) where {N<:Real}
get_b(::ResetMap{N}) where {N<:Real}
constraints_list(::ResetMap{N}) where {N<:Real}
constraints_list(::ResetMap{N, S}) where {N<:Real, S<:AbstractHyperrectangle}
```

Inherited from [`LazySet`](@ref):
* [`isbounded`](@ref isbounded(::LinearMap))
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
