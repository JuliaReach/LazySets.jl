```@meta
CurrentModule = LazySets
```

# [Ellipsoid](@id def_Ellipsoid)

```@docs
Ellipsoid
center(::Ellipsoid)
shape_matrix(::Ellipsoid)
ρ(::AbstractVector, ::Ellipsoid)
σ(::AbstractVector, ::Ellipsoid)
∈(::AbstractVector, ::Ellipsoid)
rand(::Type{Ellipsoid})
translate(::Ellipsoid, ::AbstractVector)
translate!(::Ellipsoid, ::AbstractVector)
```
Inherited from [`ConvexSet`](@ref):
* [`norm`](@ref norm(::ConvexSet, ::Real))
* [`radius`](@ref radius(::ConvexSet, ::Real))
* [`diameter`](@ref diameter(::ConvexSet, ::Real))
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
