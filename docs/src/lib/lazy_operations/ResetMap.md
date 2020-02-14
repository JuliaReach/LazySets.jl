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
get_A(::ResetMap{N}) where {N<:Real}
get_b(::ResetMap{N}) where {N<:Real}
constraints_list(::ResetMap{N}) where {N<:Real}
constraints_list(::ResetMap{N, S}) where {N<:Real, S<:AbstractHyperrectangle}
```
Inherited from [`AbstractAffineMap`](@ref):
* [`isempty`](@ref isempty(::AbstractAffineMap))
* [`isbounded`](@ref isbounded(::AbstractAffineMap))
* [`∈`](@ref ∈(::AbstractVector{N}, ::AbstractAffineMap{N}) where {N<:Real})
* [`vertices_list`](@ref vertices_list(::AbstractAffineMap{N}) where {N<:Real})
* [`linear_map`](@ref linear_map(::AbstractMatrix{N}, ::AbstractAffineMap{N}) where {N<:Real})

Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
