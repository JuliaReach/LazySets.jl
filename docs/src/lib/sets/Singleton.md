```@meta
CurrentModule = LazySets.SingletonModule
```

# [Singleton](@id def_Singleton)

```@docs
Singleton
rand(::Type{Singleton})
element(::Singleton)
translate!(::Singleton, ::AbstractVector)
rectify(S::Singleton)
project(::Singleton, ::AbstractVector{Int})
permute(::Singleton, ::AbstractVector{Int})
singleton_list(::Singleton)
linear_map(::AbstractMatrix, ::Singleton)
```

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`translate`](@ref translate(::LazySet, ::AbstractVector))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N})

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope))

Inherited from [`AbstractZonotope`](@ref):
* [`order`](@ref order(::AbstractZonotope))
* [`togrep`](@ref togrep(::AbstractZonotope))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real))
* [`constraints_list`](@ref constraints_list(::AbstractHyperrectangle{N}) where {N})

Inherited from [`AbstractSingleton`](@ref):
* [`σ`](@ref σ(::AbstractVector, ::AbstractSingleton))
* [`ρ`](@ref σ(::AbstractVector, ::AbstractSingleton))
* [`∈`](@ref ∈(::AbstractVector, ::AbstractSingleton))
* [`center`](@ref center(::AbstractSingleton))
* [`vertices`](@ref vertices(::AbstractSingleton{N}) where {N})
* [`vertices_list`](@ref vertices_list(::AbstractSingleton))
* [`radius_hyperrectangle`](@ref radius_hyperrectangle(::AbstractSingleton{N}) where {N})
* [`radius_hyperrectangle`](@ref radius_hyperrectangle(::AbstractSingleton{N}, ::Int) where {N})
* [`generators`](@ref generators(::AbstractSingleton{N}) where {N})
* [`genmat`](@ref genmat(::AbstractSingleton{N}) where {N})
* [`ngens`](@ref ngens(::AbstractSingleton))
* [`high`](@ref high(::AbstractSingleton))
* [`low`](@ref low(::AbstractSingleton))
* [`reflect`](@ref reflect(::AbstractSingleton))
