```@meta
CurrentModule = LazySets
```

# [Ellipsoid](@id def_Ellipsoid)

```@docs
Ellipsoid
ρ(::AbstractVector{N}, ::Ellipsoid{N}) where {N<:AbstractFloat}
σ(::AbstractVector{N}, ::Ellipsoid{N}) where {N<:AbstractFloat}
∈(::AbstractVector{N}, ::Ellipsoid{N}) where {N<:AbstractFloat}
rand(::Type{Ellipsoid})
center(::Ellipsoid{N}) where {N<:AbstractFloat}
translate(::Ellipsoid{N}, ::AbstractVector{N}) where {N<:AbstractFloat}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractCentrallySymmetric`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetric))
* [`isbounded`](@ref isbounded(::AbstractCentrallySymmetric))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetric))
* [`isuniversal`](@ref isuniversal(::AbstractCentrallySymmetric{N}, ::Bool=false) where {N<:Real})
* [`an_element`](@ref an_element(::AbstractCentrallySymmetric{N}) where {N<:Real})
