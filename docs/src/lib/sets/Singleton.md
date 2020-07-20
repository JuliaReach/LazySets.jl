```@meta
CurrentModule = LazySets
```

# [Singleton](@id def_Singleton)

```@docs
Singleton
rand(::Type{Singleton})
element(::Singleton{N}) where {N<:Real}
element(::Singleton{N}, ::Int) where {N<:Real}
translate(::Singleton{N}, ::AbstractVector{N}) where {N<:Real}
rectify(S::Singleton)
```
Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`singleton_list`](@ref singleton_list(::LazySet))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N<:Real})

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))

Inherited from [`AbstractZonotope`](@ref):
* [`ngens`](@ref ngens(::AbstractZonotope))
* [`order`](@ref order(::AbstractZonotope))
* [`togrep`](@ref togrep(::AbstractZonotope))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real))
* [`high`](@ref high(::AbstractHyperrectangle{N}) where {N<:Real})
* [`low`](@ref low(::AbstractHyperrectangle{N}) where {N<:Real})
* [`constraints_list`](@ref constraints_list(::AbstractHyperrectangle{N}) where {N<:Real})

Inherited from [`AbstractSingleton`](@ref):
* [`σ`](@ref σ(::AbstractVector{N}, ::AbstractSingleton{N}) where {N<:Real})
* [`∈`](@ref ∈(::AbstractVector{N}, ::AbstractSingleton{N}) where {N<:Real})
* [`an_element`](@ref an_element(::AbstractSingleton{N}) where {N<:Real})
* [`center`](@ref center(::AbstractSingleton{N}) where {N<:Real})
* [`vertices_list`](@ref vertices_list(::AbstractSingleton{N}) where {N<:Real})
* [`radius_hyperrectangle`](@ref radius_hyperrectangle(::AbstractSingleton{N}) where {N<:Real})
* [`radius_hyperrectangle`](@ref radius_hyperrectangle(::AbstractSingleton{N}, ::Int) where {N<:Real})
* [`linear_map`](@ref linear_map(::AbstractMatrix{N}, ::AbstractSingleton{N}) where {N<:Real})
* [`generators`](@ref generators(::AbstractSingleton{N}) where {N<:Real})
* [`genmat`](@ref genmat(::AbstractSingleton{N}) where {N<:Real})
