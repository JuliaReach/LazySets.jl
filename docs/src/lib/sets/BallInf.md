```@meta
CurrentModule = LazySets
```

# [Infinity-norm ball (BallInf)](@id def_BallInf)

```@docs
BallInf
center(::BallInf)
radius_ball(::BallInf)
radius(::BallInf, ::Real=Inf)
radius_hyperrectangle(::BallInf)
radius_hyperrectangle(::BallInf, ::Int)
isflat(::BallInf)
rand(::Type{BallInf})
σ(::AbstractVector, ::BallInf)
ρ(::AbstractVector, ::BallInf)
translate(::BallInf, ::AbstractVector)
translate!(::BallInf, ::AbstractVector)
ngens(::BallInf)
volume(::BallInf)
reflect(::BallInf)
```
Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`area`](@ref area(::LazySet))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N})

Inherited from [`AbstractZonotope`](@ref):
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractZonotope))
* [`order`](@ref order(::AbstractZonotope))
* [`togrep`](@ref togrep(::AbstractZonotope))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`∈`](@ref ∈(::AbstractVector, ::AbstractHyperrectangle))
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`vertices_list`](@ref vertices_list(::AbstractHyperrectangle))
* [`high`](@ref high(::AbstractHyperrectangle))
* [`high`](@ref high(::AbstractHyperrectangle, ::Int))
* [`low`](@ref low(::AbstractHyperrectangle))
* [`low`](@ref low(::AbstractHyperrectangle, ::Int))
* [`extrema`](@ref extrema(::AbstractHyperrectangle))
* [`extrema`](@ref extrema(::AbstractHyperrectangle, ::Int))
* [`generators`](@ref generators(::AbstractHyperrectangle))
* [`genmat`](@ref genmat(::AbstractHyperrectangle))
* [`constraints_list`](@ref constraints_list(::AbstractHyperrectangle{N}) where {N})
* [`rectify`](@ref rectify(::AbstractHyperrectangle))
