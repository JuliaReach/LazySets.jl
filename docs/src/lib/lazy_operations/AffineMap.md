```@meta
CurrentModule = LazySets
```

# [Affine map (AffineMap)](@id def_AffineMap)

```@docs
AffineMap
```

Inherited from [`AbstractAffineMap`](@ref):
* [`dim`](@ref dim(::AbstractAffineMap))
* [`σ`](@ref σ(::AbstractVector{N}, ::AbstractAffineMap{N}) where {N<:Real})
* [`ρ`](@ref ρ(::AbstractVector{N}, ::AbstractAffineMap{N}) where {N<:Real})
* [`an_element`](@ref an_element(::AbstractAffineMap))
* [`isempty`](@ref isempty(::AbstractAffineMap))
* [`isbounded`](@ref isbounded(::AbstractAffineMap))
* [`∈`](@ref ∈(::AbstractVector{N}, ::AbstractAffineMap{N}) where {N<:Real})
* [`vertices_list`](@ref vertices_list(::AbstractAffineMap{N}) where {N<:Real})
* [`constraints_list`](@ref constraints_list(::AbstractAffineMap{N}) where {N<:Real})
* [`linear_map`](@ref linear_map(::AbstractMatrix{N}, ::AbstractAffineMap{N}) where {N<:Real})

Inherited from [`LazySet`](@ref):
* [`singleton_list`](@ref singleton_list(::LazySet))
