```@meta
CurrentModule = LazySets
```

# [Origin (ZeroSet)](@id def_ZeroSet)

```@docs
ZeroSet
dim(::ZeroSet)
σ(::AbstractVector{N}, ::ZeroSet{N}) where {N<:Real}
ρ(::AbstractVector{N}, ::ZeroSet{N}) where {N<:Real}
∈(::AbstractVector{N}, ::ZeroSet{N}) where {N<:Real}
rand(::Type{ZeroSet})
element(::ZeroSet{N}) where {N<:Real}
element(::ZeroSet{N}, ::Int) where {N<:Real}
linear_map(::AbstractMatrix{N}, ::ZeroSet{N}) where {N<:Real}
translate(::ZeroSet{N}, ::AbstractVector{N}) where {N<:Real}
center(::ZeroSet{N}, ::Int) where {N<:Real}
rectify(Z::ZeroSet)
```
Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`singleton_list`](@ref singleton_list(::LazySet))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N})

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
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
* [`radius_hyperrectangle`](@ref radius_hyperrectangle(::AbstractSingleton{N}) where {N})
* [`radius_hyperrectangle`](@ref radius_hyperrectangle(::AbstractSingleton{N}, ::Int) where {N})
* [`vertices`](@ref vertices(::AbstractSingleton{N}) where {N})
* [`vertices_list`](@ref vertices_list(::AbstractSingleton))
* [`center`](@ref center(::AbstractSingleton))
* [`generators`](@ref generators(::AbstractSingleton{N}) where {N})
* [`genmat`](@ref genmat(::AbstractSingleton{N}) where {N})
