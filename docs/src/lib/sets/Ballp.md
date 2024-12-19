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
```@meta
CurrentModule = LazySets
```
* [`ball_norm`](@ref ball_norm(::LazySet))
* [`radius_ball`](@ref radius_ball(::LazySet))

Inherited from [`LazySet`](@ref):
* [`concretize`](@ref concretize(::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`eltype`](@ref eltype(::LazySet))
* [`high`](@ref high(::LazySet))
* [`low`](@ref low(::LazySet))
* [`ispolyhedral`](@ref ispolyhedral(::LazySet))
* [`isoperation`](@ref isoperation(::LazySet))
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`rectify`](@ref rectify(::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`sample`](@ref sample(::LazySet, ::Int))
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
* [`≈`](@ref ≈(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))

Inherited from [`AbstractCentrallySymmetric`](@ref):
* [`an_element`](@ref an_element(::AbstractCentrallySymmetric))
* [`center`](@ref center(::AbstractCentrallySymmetric, ::Int))
* [`dim`](@ref dim(::AbstractCentrallySymmetric))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetric))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetric, ::Int))
* [`isbounded`](@ref isbounded(::AbstractCentrallySymmetric))
* [`isboundedtype`](@ref isboundedtype(::Type{AbstractCentrallySymmetric}))
* [`isconvextype`](@ref isconvextype(::Type{AbstractCentrallySymmetric}))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetric))
* [`isuniversal`](@ref isuniversal(::AbstractCentrallySymmetric, ::Bool=false))

Inherited from [`AbstractBallp`](@ref):
* [`low`](@ref low(::AbstractBallp, ::Int))
* [`high`](@ref high(::AbstractBallp, ::Int))
* [`∈`](@ref ∈(::AbstractVector, ::AbstractBallp))
* [`ρ`](@ref ρ(::AbstractVector, ::AbstractBallp))
* [`σ`](@ref σ(::AbstractVector, ::AbstractBallp))
```@meta
CurrentModule = LazySets.API
```
* [`minkowski_sum`](@ref minkowski_sum(::LazySet, ::LazySet))
