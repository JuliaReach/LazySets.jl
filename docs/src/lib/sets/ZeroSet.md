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
```
Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N<:Real})
* [`singleton_list`](@ref singleton_list(::AbstractPolytope{N}) where {N<:Real})

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
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
* [`radius_hyperrectangle`](@ref radius_hyperrectangle(::AbstractSingleton{N}) where {N<:Real})
* [`radius_hyperrectangle`](@ref radius_hyperrectangle(::AbstractSingleton{N}, ::Int) where {N<:Real})
* [`vertices_list`](@ref vertices_list(::AbstractSingleton{N}) where {N<:Real})
* [`center`](@ref center(::AbstractSingleton{N}) where {N<:Real})
* [`an_element`](@ref an_element(::AbstractSingleton{N}) where {N<:Real})
* [`generators`](@ref generators(::AbstractSingleton{N}) where {N<:Real})
* [`genmat`](@ref genmat(::AbstractSingleton{N}) where {N<:Real})
