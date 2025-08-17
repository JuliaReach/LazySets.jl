```@meta
CurrentModule = LazySets.BallpModule
```

# [p-norm ball (Ballp)](@id def_Ballp)

```@docs
Ballp
```

## Operations

```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
rand(::Type{LazySet})
```
```@meta
CurrentModule = LazySets.BallpModule
```
```@docs
rand(::Type{Ballp})
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
reflect(::LazySet)
```
```@meta
CurrentModule = LazySets.BallpModule
```
```@docs
reflect(::Ballp)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
translate!(::LazySet, ::AbstractVector)
```
```@meta
CurrentModule = LazySets.BallpModule
```
```@docs
translate!(::Ballp, ::AbstractVector)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:
* [`center`](@ref center(::LazySet))
* [`isoperationtype`](@ref isoperationtype(::Type{LazySet}))
* [`project`](@ref project(::LazySet, ::AbstractVector{Int}))
* [`scale`](@ref scale(::Real, ::LazySet))
```@meta
CurrentModule = LazySets
```
* [`ball_norm`](@ref ball_norm(::LazySet))
* [`radius_ball`](@ref radius_ball(::LazySet))

Inherited from [`LazySet`](@ref):
* [`chebyshev_center_radius`](@ref chebyshev_center_radius(::LazySet))
* [`concretize`](@ref concretize(::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet))
* `copy(::Type{LazySet})`
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`eltype`](@ref eltype(::LazySet))
* [`high`](@ref high(::LazySet))
* [`isoperation`](@ref isoperation(::LazySet))
* [`ispolyhedral`](@ref ispolyhedral(::LazySet))
* [`ispolyhedraltype`](@ref ispolyhedraltype(::Type{LazySet}))
* [`ispolytopic`](@ref ispolytopic(::LazySet))
* [`low`](@ref low(::LazySet))
* [`norm`](@ref norm(::LazySet, ::Real))
* [`polyhedron`](@ref polyhedron(::LazySet))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`rationalize`](@ref rationalize(::LazySet))
* [`rectify`](@ref rectify(::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`sample`](@ref sample(::LazySet, ::Int))
* [`tosimplehrep`](@ref tosimplehrep(::LazySet))
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
* [`triangulate`](@ref triangulate(::LazySet))
* [`triangulate_faces`](@ref triangulate_faces(::LazySet))
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

Inherited from [`AbstractBallp`](@ref):
* [`low`](@ref low(::AbstractBallp, ::Int))
* [`high`](@ref high(::AbstractBallp, ::Int))
* [`∈`](@ref ∈(::AbstractVector, ::AbstractBallp))
* [`ρ`](@ref ρ(::AbstractVector, ::AbstractBallp))
* [`σ`](@ref σ(::AbstractVector, ::AbstractBallp))
* [`minkowski_sum`](@ref minkowski_sum(::AbstractBallp, ::AbstractBallp))
