```@meta
CurrentModule = LazySets.Ball2Module
```

# [Euclidean-norm ball (Ball2)](@id def_Ball2)

```@docs
Ball2
```

## Operations

```@docs
chebyshev_center_radius(::Ball2)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
rand(::Type{LazySet})
```
```@meta
CurrentModule = LazySets.Ball2Module
```
```@docs
rand(::Type{Ball2})
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
reflect(::LazySet)
```
```@meta
CurrentModule = LazySets.Ball2Module
```
```@docs
reflect(::Ball2)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
volume(::LazySet)
```
```@meta
CurrentModule = LazySets.Ball2Module
```
```@docs
volume(::Ball2)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
∈(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets.Ball2Module
```
```@docs
∈(::AbstractVector, ::Ball2)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
sample(::LazySet, ::Int)
```
```@meta
CurrentModule = LazySets.Ball2Module
```
```@docs
sample(::Ball2, ::Int)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
ρ(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets.Ball2Module
```
```@docs
ρ(::AbstractVector, ::Ball2)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
σ(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets.Ball2Module
```
```@docs
σ(::AbstractVector, ::Ball2)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
translate!(::LazySet, ::AbstractVector)
```
```@meta
CurrentModule = LazySets.Ball2Module
```
```@docs
translate!(::Ball2, ::AbstractVector)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
isdisjoint(::LazySet, ::LazySet, ::Bool=false)
```
```@meta
CurrentModule = LazySets.Ball2Module
```
```@docs
isdisjoint(::Ball2, ::Ball2, ::Bool=false)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
⊆(::LazySet, ::LazySet, ::Bool=false)
```
```@meta
CurrentModule = LazySets.Ball2Module
```
```@docs
⊆(::Ball2, ::Ball2, ::Bool=false)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:
* [`area`](@ref area(::LazySet))
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
* [`concretize`](@ref concretize(::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`eltype`](@ref eltype(::LazySet))
* [`high`](@ref high(::LazySet))
* [`ispolyhedral`](@ref ispolyhedral(::LazySet))
* [`low`](@ref low(::LazySet))
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`rectify`](@ref rectify(::LazySet))
* [`surface`](@ref surface(::LazySet))
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`cartesian_product`](@ref cartesian_product(::LazySet, ::LazySet))
* [`≈`](@ref ≈(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))

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
* [`isuniversal`](@ref isuniversal(::AbstractCentrallySymmetric, ::Bool=false))

Inherited from [`AbstractBallp`](@ref):
* [`high`](@ref high(::AbstractBallp, ::Int))
* [`low`](@ref low(::AbstractBallp, ::Int))
* `minkowski_sum`

```@meta
CurrentModule = LazySets.Ball2Module
```

```@docs
_sample_unit_nball_muller!
```
