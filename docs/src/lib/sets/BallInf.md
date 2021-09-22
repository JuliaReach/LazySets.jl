```@meta
CurrentModule = LazySets
```

# [Infinity-norm ball (BallInf)](@id def_BallInf)

```@docs
BallInf
center(::BallInf)
radius(::BallInf, ::Real=Inf)
radius_hyperrectangle(::BallInf)
radius_hyperrectangle(::BallInf, ::Int)
isflat(::BallInf)
rand(::Type{BallInf})
σ(::AbstractVector, ::BallInf)
ρ(::AbstractVector, ::BallInf)
translate(::BallInf, ::AbstractVector)
translate!(::BallInf, ::AbstractVector)
```
Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`singleton_list`](@ref singleton_list(::LazySet))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N})

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope))

Inherited from [`AbstractZonotope`](@ref):
* [`ngens`](@ref ngens(::AbstractZonotope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractZonotope))
* [`order`](@ref order(::AbstractZonotope))
* [`togrep`](@ref togrep(::AbstractZonotope))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`∈`](@ref ∈(::AbstractVector, ::AbstractHyperrectangle))
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`vertices_list`](@ref vertices_list(::AbstractHyperrectangle))
* [`high`](@ref high(::AbstractHyperrectangle))
* [`low`](@ref low(::AbstractHyperrectangle))
* [`generators`](@ref generators(::AbstractHyperrectangle))
* [`genmat`](@ref genmat(::AbstractHyperrectangle))
* [`constraints_list`](@ref constraints_list(::AbstractHyperrectangle{N}) where {N})
