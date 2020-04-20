```@meta
CurrentModule = LazySets
```

# [p-norm ball (Ballp)](@id def_Ballp)

```@docs
Ballp
σ(::AbstractVector{N}, ::Ballp{N}) where {N<:AbstractFloat}
∈(::AbstractVector{N}, ::Ballp{N}) where {N<:AbstractFloat}
center(::Ballp{N}) where {N<:AbstractFloat}
rand(::Type{Ballp})
translate(::Ballp{N}, ::AbstractVector{N}) where {N<:AbstractFloat}
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
