```@meta
CurrentModule = LazySets.EllipsoidModule
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
translate!(::Ellipsoid, ::AbstractVector)
linear_map(::AbstractMatrix, ::Ellipsoid)
```

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`rectify`](@ref rectify(::LazySet))
* [`low`](@ref low(::LazySet))
* [`high`](@ref high(::LazySet))
* [`translate`](@ref translate(::LazySet, ::AbstractVector))

Inherited from [`AbstractCentrallySymmetric`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetric))
* [`isbounded`](@ref isbounded(::AbstractCentrallySymmetric))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetric))
* [`isuniversal`](@ref isuniversal(::AbstractCentrallySymmetric{N}, ::Bool=false) where {N})
* [`an_element`](@ref an_element(::AbstractCentrallySymmetric))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetric))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetric, ::Int))
