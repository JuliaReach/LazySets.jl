```@meta
CurrentModule = LazySets
```

# [Singleton](@id def_Singleton)

```@docs
Singleton
rand(::Type{Singleton})
element(::Singleton)
element(::Singleton, ::Int)
translate(::Singleton, ::AbstractVector)
rectify(S::Singleton)
project(::Singleton, ::AbstractVector{Int})
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
* [`order`](@ref order(::AbstractZonotope))
* [`togrep`](@ref togrep(::AbstractZonotope))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real))
* [`high`](@ref high(::AbstractHyperrectangle))
* [`low`](@ref low(::AbstractHyperrectangle))
* [`constraints_list`](@ref constraints_list(::AbstractHyperrectangle{N}) where {N})

Inherited from [`AbstractSingleton`](@ref):
* [`σ`](@ref σ(::AbstractVector, ::AbstractSingleton))
* [`∈`](@ref ∈(::AbstractVector, ::AbstractSingleton))
* [`center`](@ref center(::AbstractSingleton))
* [`vertices`](@ref vertices(::AbstractSingleton{N}) where {N})
* [`vertices_list`](@ref vertices_list(::AbstractSingleton))
* [`radius_hyperrectangle`](@ref radius_hyperrectangle(::AbstractSingleton{N}) where {N})
* [`radius_hyperrectangle`](@ref radius_hyperrectangle(::AbstractSingleton{N}, ::Int) where {N})
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractSingleton))
* [`generators`](@ref generators(::AbstractSingleton{N}) where {N})
* [`genmat`](@ref genmat(::AbstractSingleton{N}) where {N})
