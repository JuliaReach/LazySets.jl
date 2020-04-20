```@meta
CurrentModule = LazySets
```

# [Infinity-norm ball (BallInf)](@id def_BallInf)

```@docs
BallInf
center(::BallInf{N}) where {N<:Real}
radius(::BallInf, ::Real=Inf)
radius_hyperrectangle(::BallInf{N}) where {N<:Real}
radius_hyperrectangle(::BallInf{N}, ::Int) where {N<:Real}
isflat(::BallInf)
rand(::Type{BallInf})
σ(::AbstractVector{N}, ::BallInf{N}) where {N<:Real}
ρ(::AbstractVector{N}, ::BallInf{N}) where {N<:Real}
translate(::BallInf{N}, ::AbstractVector{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N<:Real})
* [`singleton_list`](@ref singleton_list(::AbstractPolytope{N}) where {N<:Real})

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope{N}) where {N<:Real})

Inherited from [`AbstractZonotope`](@ref):
* [`ngens`](@ref ngens(::AbstractZonotope))
* [`linear_map`](@ref linear_map(::AbstractMatrix{N}, ::AbstractZonotope{N}) where {N<:Real})
* [`order`](@ref order(::AbstractZonotope))
* [`togrep`](@ref togrep(::AbstractZonotope))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`∈`](@ref ∈(::AbstractVector{N}, ::AbstractHyperrectangle{N}) where {N<:Real})
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`vertices_list`](@ref vertices_list(::AbstractHyperrectangle{N}) where {N<:Real})
* [`high`](@ref high(::AbstractHyperrectangle{N}) where {N<:Real})
* [`low`](@ref low(::AbstractHyperrectangle{N}) where {N<:Real})
* [`generators`](@ref generators(::AbstractHyperrectangle))
* [`genmat`](@ref genmat(::AbstractHyperrectangle))
* [`constraints_list`](@ref constraints_list(::AbstractHyperrectangle{N}) where {N<:Real})
