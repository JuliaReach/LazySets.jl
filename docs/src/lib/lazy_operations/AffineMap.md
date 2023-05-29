```@meta
CurrentModule = LazySets
```

# [Affine map (AffineMap)](@id def_AffineMap)

```@docs
AffineMap
```

Inherited from [`AbstractAffineMap`](@ref):
* [`dim`](@ref dim(::AbstractAffineMap))
* [`σ`](@ref σ(::AbstractVector, ::AbstractAffineMap))
* [`ρ`](@ref ρ(::AbstractVector, ::AbstractAffineMap))
* [`an_element`](@ref an_element(::AbstractAffineMap))
* [`isempty`](@ref isempty(::AbstractAffineMap))
* [`isbounded`](@ref isbounded(::AbstractAffineMap))
* [`∈`](@ref ∈(::AbstractVector, ::AbstractAffineMap))
* [`vertices_list`](@ref vertices_list(::AbstractAffineMap))
* [`constraints_list`](@ref constraints_list(::AbstractAffineMap))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractAffineMap))

Inherited from [`LazySet`](@ref):
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`reflect`](@ref reflect(::LazySet))
