```@meta
CurrentModule = LazySets
```

# [Symmetric interval hull (SymmetricIntervalHull)](@id def_SymmetricIntervalHull)

```@docs
SymmetricIntervalHull
dim(::SymmetricIntervalHull)
σ(::AbstractVector{N}, ::SymmetricIntervalHull{N}) where {N<:Real}
center(::SymmetricIntervalHull{N}, ::Int) where {N<:Real}
center(::SymmetricIntervalHull{N}) where {N<:Real}
radius_hyperrectangle(::SymmetricIntervalHull{N}) where {N<:Real}
radius_hyperrectangle(::SymmetricIntervalHull{N}, ::Int) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`singleton_list`](@ref singleton_list(::AbstractPolytope{N}) where {N<:Real})

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope{N}) where {N<:Real})

Inherited from [`AbstractZonotope`](@ref):
* [`ngens`](@ref ngens(::AbstractZonotope))
* [`linear_map`](@ref linear_map(::AbstractMatrix{N}, ::AbstractZonotope{N}) where {N<:Real})
* [`order`](@ref order(::AbstractZonotope))
* [`togrep`](@ref togrep(::AbstractZonotope))
* [`translate`](@ref translate(::AbstractZonotope{N}, ::AbstractVector{N}) where {N<:Real})

Inherited from [`AbstractHyperrectangle`](@ref):
* [`ρ`](@ref ρ(::AbstractVector{N}, ::AbstractHyperrectangle{N}) where {N<:Real})
* [`σ`](@ref σ(::AbstractVector{N}, ::AbstractHyperrectangle{N}) where {N<:Real})
* [`∈`](@ref ∈(::AbstractVector{N}, ::AbstractHyperrectangle{N}) where {N<:Real})
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real))
* [`constraints_list`](@ref constraints_list(::AbstractHyperrectangle{N}) where {N<:Real})
* [`vertices_list`](@ref vertices_list(::AbstractHyperrectangle{N}) where {N<:Real})
* [`high`](@ref high(::AbstractHyperrectangle{N}) where {N<:Real})
* [`low`](@ref low(::AbstractHyperrectangle{N}) where {N<:Real})
* [`generators`](@ref generators(::AbstractHyperrectangle))
* [`genmat`](@ref genmat(::AbstractHyperrectangle))
