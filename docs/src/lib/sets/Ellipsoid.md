```@meta
CurrentModule = LazySets.EllipsoidModule
```

# [Ellipsoid](@id def_Ellipsoid)

```@docs
Ellipsoid
```

## Operations

```@docs
center(::Ellipsoid)
rand(::Type{Ellipsoid})
shape_matrix(::Ellipsoid)
∈(::AbstractVector, ::Ellipsoid)
ρ(::AbstractVector, ::Ellipsoid)
σ(::AbstractVector, ::Ellipsoid)
linear_map(::AbstractMatrix, ::Ellipsoid)
translate!(::Ellipsoid, ::AbstractVector)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:
* [`isoperationtype`](@ref isoperationtype(::Type{LazySet}))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`concretize`](@ref concretize(::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`eltype`](@ref eltype(::LazySet))
* [`high`](@ref high(::LazySet))
* [`is_polyhedral`](@ref is_polyhedral(::LazySet))
* [`isoperation`](@ref isoperation(::LazySet))
* [`low`](@ref low(::LazySet))
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`sample`](@ref sample(::LazySet, ::Int))
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
* [`≈`](@ref ≈(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))

Inherited from [`ConvexSet`](@ref):
* [`high`](@ref high(::ConvexSet, ::Int))
* [`low`](@ref low(::ConvexSet, ::Int))

Inherited from [`AbstractCentrallySymmetric`](@ref):
* [`an_element`](@ref an_element(::AbstractCentrallySymmetric))
* [`center`](@ref center(::AbstractCentrallySymmetric, ::Int))
* [`dim`](@ref dim(::AbstractCentrallySymmetric))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetric))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetric, ::Int))
* [`isbounded`](@ref isbounded(::AbstractCentrallySymmetric))
```@meta
CurrentModule = LazySets.API
```
* [`isboundedtype`](@ref isboundedtype(::Type{LazySet}))
* [`isconvextype`](@ref isconvextype(::Type{LazySet}))
```@meta
CurrentModule = LazySets
```
* [`isempty`](@ref isempty(::AbstractCentrallySymmetric))
* [`isuniversal`](@ref isuniversal(::AbstractCentrallySymmetric{N}, ::Bool=false) where {N})
