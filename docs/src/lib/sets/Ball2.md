```@meta
CurrentModule = LazySets
```

# [Euclidean-norm ball (Ball2)](@id def_Ball2)

```@docs
Ball2
σ(::AbstractVector{N}, ::Ball2{N}) where {N<:AbstractFloat}
∈(::AbstractVector{N}, ::Ball2{N}) where {N<:AbstractFloat}
center(::Ball2{N}) where {N<:AbstractFloat}
rand(::Type{Ball2})
sample(::Ball2{N}, ::Int) where {N<:AbstractFloat}
center(::Ball2{N}, ::Int) where {N<:Real}
translate(::Ball2{N}, ::AbstractVector{N}) where {N<:AbstractFloat}
chebyshev_center(::Ball2{N}) where {N<:AbstractFloat}
volume(::Ball2{N}) where {N<:AbstractFloat}
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
