```@meta
CurrentModule = LazySets
```

# [Symmetric interval hull (SymmetricIntervalHull)](@id def_SymmetricIntervalHull)

```@docs
SymmetricIntervalHull
⊡
dim(::SymmetricIntervalHull)
σ(::AbstractVector, ::SymmetricIntervalHull)
center(::SymmetricIntervalHull{N}, ::Int) where {N}
center(::SymmetricIntervalHull{N}) where {N}
radius_hyperrectangle(::SymmetricIntervalHull)
radius_hyperrectangle(::SymmetricIntervalHull, ::Int)
```
Inherited from [`ConvexSet`](@ref):
* [`diameter`](@ref diameter(::ConvexSet, ::Real))
* [`singleton_list`](@ref singleton_list(::ConvexSet))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope))

Inherited from [`AbstractZonotope`](@ref):
* [`ngens`](@ref ngens(::AbstractZonotope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractZonotope))
* [`order`](@ref order(::AbstractZonotope))
* [`togrep`](@ref togrep(::AbstractZonotope))
* [`translate`](@ref translate(::AbstractZonotope, ::AbstractVector))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`ρ`](@ref ρ(::AbstractVector, ::AbstractHyperrectangle))
* [`σ`](@ref σ(::AbstractVector, ::AbstractHyperrectangle))
* [`∈`](@ref ∈(::AbstractVector, ::AbstractHyperrectangle))
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real))
* [`constraints_list`](@ref constraints_list(::AbstractHyperrectangle))
* [`vertices_list`](@ref vertices_list(::AbstractHyperrectangle))
* [`high`](@ref high(::AbstractHyperrectangle))
* [`low`](@ref low(::AbstractHyperrectangle))
* [`generators`](@ref generators(::AbstractHyperrectangle))
* [`genmat`](@ref genmat(::AbstractHyperrectangle))
