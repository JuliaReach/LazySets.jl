```@meta
CurrentModule = LazySets.Ball2Module
```

# [Euclidean-norm ball (Ball2)](@id def_Ball2)

```@docs
Ball2
```

## Operations

```@docs
center(::Ball2)
chebyshev_center_radius(::Ball2)
rand(::Type{Ball2})
reflect(::Ball2)
volume(::Ball2)
∈(::AbstractVector, ::Ball2)
sample(::Ball2{N}, ::Int) where {N}
ρ(::AbstractVector, ::Ball2)
σ(::AbstractVector, ::Ball2)
translate!(::Ball2, ::AbstractVector)
isdisjoint(::Ball2, ::Ball2, ::Bool=false)
⊆(::Ball2, ::Ball2, ::Bool=false)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:
* [`area`](@ref area(::LazySet))
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
* [`isuniversal`](@ref isuniversal(::AbstractCentrallySymmetric{N}, ::Bool=false) where {N})

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
