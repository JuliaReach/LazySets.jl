```@meta
CurrentModule = LazySets
```

# [Ellipsoid](@id def_Ellipsoid)

```@docs
Ellipsoid
ρ(::AbstractVector, ::Ellipsoid)
σ(::AbstractVector, ::Ellipsoid)
∈(::AbstractVector, ::Ellipsoid)
rand(::Type{Ellipsoid})
center(::Ellipsoid)
shape_matrix(::Ellipsoid)
translate(::Ellipsoid, ::AbstractVector)
translate!(::Ellipsoid, ::AbstractVector)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractCentrallySymmetric`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetric))
* [`isbounded`](@ref isbounded(::AbstractCentrallySymmetric))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetric))
* [`isuniversal`](@ref isuniversal(::AbstractCentrallySymmetric{N}, ::Bool=false) where {N})
* [`an_element`](@ref an_element(::AbstractCentrallySymmetric))
