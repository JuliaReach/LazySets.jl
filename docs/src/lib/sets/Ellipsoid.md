```@meta
CurrentModule = LazySets.EllipsoidModule
```

# [Ellipsoid](@id def_Ellipsoid)

```@docs
Ellipsoid
```

## Operations

```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
rand(::Type{LazySet})
```
```@meta
CurrentModule = LazySets.EllipsoidModule
```
```@docs
rand(::Type{Ellipsoid})
shape_matrix(::Ellipsoid)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
∈(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets.EllipsoidModule
```
```@docs
∈(::AbstractVector, ::Ellipsoid)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
ρ(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets.EllipsoidModule
```
```@docs
ρ(::AbstractVector, ::Ellipsoid)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
σ(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets.EllipsoidModule
```
```@docs
σ(::AbstractVector, ::Ellipsoid)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
linear_map(::AbstractMatrix, ::LazySet)
```
```@meta
CurrentModule = LazySets.EllipsoidModule
```
```@docs
linear_map(::AbstractMatrix, ::Ellipsoid)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:
* [`center`](@ref center(::LazySet))
* [`isoperationtype`](@ref isoperationtype(::Type{LazySet}))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`translate!`](@ref translate!(::LazySet, ::AbstractVector))

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`chebyshev_center_radius`](@ref chebyshev_center_radius(::LazySet))
* [`concretize`](@ref concretize(::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet))
* `copy(::Type{LazySet})`
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`eltype`](@ref eltype(::LazySet))
* [`high`](@ref high(::LazySet))
* [`high`](@ref high(::LazySet, ::Int))
* [`ispolyhedral`](@ref ispolyhedral(::LazySet))
* [`isoperation`](@ref isoperation(::LazySet))
* [`low`](@ref low(::LazySet))
* [`low`](@ref low(::LazySet, ::Int))
* [`norm`](@ref norm(::LazySet, ::Real))
* [`polyhedron`](@ref polyhedron(::LazySet))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`rationalize`](@ref rationalize(::LazySet))
* [`triangulate`](@ref triangulate(::LazySet))
* [`triangulate_faces`](@ref triangulate_faces(::LazySet))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`project`](@ref project(::LazySet, ::AbstractVector{Int}))
* [`sample`](@ref sample(::LazySet, ::Int))
* [`tosimplehrep`](@ref tosimplehrep(::LazySet))
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
* [`convex_hull`](@ref convex_hull(::LazySet, ::LazySet))
* [`isapprox`](@ref isapprox(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))

Inherited from [`ConvexSet`](@ref):
* [`isconvextype`](@ref isconvextype(::Type{ConvexSet}))
* [`linear_combination`](@ref linear_combination(::ConvexSet, ::ConvexSet))

Inherited from [`AbstractCentrallySymmetric`](@ref):
* [`an_element`](@ref an_element(::AbstractCentrallySymmetric))
* [`center`](@ref center(::AbstractCentrallySymmetric, ::Int))
* [`dim`](@ref dim(::AbstractCentrallySymmetric))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetric))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetric, ::Int))
* [`isbounded`](@ref isbounded(::AbstractCentrallySymmetric))
* [`isboundedtype`](@ref isboundedtype(::Type{AbstractCentrallySymmetric}))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetric))
* [`isuniversal`](@ref isuniversal(::AbstractCentrallySymmetric, ::Bool=false))
