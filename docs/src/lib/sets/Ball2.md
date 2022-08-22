```@meta
CurrentModule = LazySets
```

# [Euclidean-norm ball (Ball2)](@id def_Ball2)

```@docs
Ball2
ρ(::AbstractVector, ::Ball2)
σ(::AbstractVector, ::Ball2)
∈(::AbstractVector, ::Ball2)
center(::Ball2)
rand(::Type{Ball2})
sample(::Ball2{N}, ::Int) where {N}
translate(::Ball2, ::AbstractVector)
translate!(::Ball2, ::AbstractVector)
chebyshev_center_radius(::Ball2)
volume(::Ball2)
```
Inherited from [`ConvexSet`](@ref):
* [`norm`](@ref norm(::ConvexSet, ::Real))
* [`radius`](@ref radius(::ConvexSet, ::Real))
* [`diameter`](@ref diameter(::ConvexSet, ::Real))
* [`area`](@ref area(::ConvexSet))
* [`rectify`](@ref rectify(::ConvexSet))
* [`low`](@ref low(::ConvexSet))
* [`low`](@ref low(::ConvexSet{N}, ::Int) where {N})
* [`high`](@ref high(::ConvexSet))
* [`high`](@ref high(::ConvexSet{N}, ::Int) where {N})

Inherited from [`AbstractCentrallySymmetric`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetric))
* [`isbounded`](@ref isbounded(::AbstractCentrallySymmetric))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetric))
* [`isuniversal`](@ref isuniversal(::AbstractCentrallySymmetric{N}, ::Bool=false) where {N})
* [`an_element`](@ref an_element(::AbstractCentrallySymmetric))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetric))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetric, ::Int))
