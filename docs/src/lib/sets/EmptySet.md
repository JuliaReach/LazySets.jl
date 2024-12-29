```@meta
CurrentModule = LazySets.EmptySetModule
```

# [Empty set (EmptySet)](@id def_EmptySet)

```@docs
EmptySet
∅
```

## Conversion

```julia
convert(::Type{EmptySet}, ::LazySet)
```

## Operations

```@docs
chebyshev_center_radius(::EmptySet; kwargs...)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
complement(::LazySet)
```
```@meta
CurrentModule = LazySets.EmptySetModule
```
```@docs
complement(::EmptySet)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
rand(::Type{<:LazySet})
```
```@meta
CurrentModule = LazySets.EmptySetModule
```
```@docs
rand(::Type{EmptySet})
plot_recipe(::EmptySet{N}, ::Any=zero(N)) where {N}
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:

* [`an_element`](@ref an_element(::LazySet))
* [`area`](@ref area(::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet))
* `copy(::EmptySet)`
* [`diameter`](@ref diameter(::LazySet, ::Real=Inf))
* [`dim`](@ref dim(::LazySet))
* [`high`](@ref high(::LazySet))
* [`high`](@ref high(::LazySet, ::Int))
* [`isbounded`](@ref isbounded(::LazySet))
* [`isboundedtype`](@ref isboundedtype(::Type{LazySet}))
* [`isconvextype`](@ref isconvextype(::Type{LazySet}))
* [`isempty`](@ref isempty(::LazySet))
* [`isoperationtype`](@ref isoperationtype(::Type{LazySet}))
* [`isuniversal`](@ref isuniversal(::LazySet, ::Bool=false))
* [`low`](@ref low(::LazySet))
* [`low`](@ref low(::LazySet, ::Int))
* [`norm`](@ref norm(::LazySet, ::Real=Inf))
* [`radius`](@ref radius(::LazySet, ::Real=Inf))
* [`rationalize`](@ref rationalize(::LazySet))
* [`rectify`](@ref rectify(::LazySet))
* [`reflect`](@ref reflect(::LazySet))
* [`surface`](@ref surface(::LazySet))
* [`vertices_list`](@ref vertices_list(::LazySet))
* [`vertices`](@ref vertices(::LazySet))
* [`volume`](@ref volume(::LazySet))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`∈`](@ref ∈(::AbstractVector, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::LazySet))
* [`permute`](@ref permute(::LazySet, ::AbstractVector{Int}))
* [`project`](@ref project(::LazySet, ::AbstractVector))
* [`sample`](@ref sample(::LazySet, ::Int=1))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`scale!`](@ref scale!(::Real, ::LazySet))
* [`ρ`](@ref ρ(::AbstractVector, ::LazySet))
* [`σ`](@ref σ(::AbstractVector, ::LazySet))
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
* [`translate!`](@ref translate!(::LazySet, ::AbstractVector))
* [`cartesian_product`](@ref cartesian_product(::LazySet, ::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet, ::LazySet))
* [`difference`](@ref difference(::LazySet, ::LazySet))
* [`distance`](@ref distance(::LazySet, ::LazySet; ::Real=2.0))
* [`intersection`](@ref intersection(::LazySet, ::LazySet))
* [`≈`](@ref ≈(::LazySet, ::LazySet))
* [`isdisjoint`](@ref isdisjoint(::LazySet, ::LazySet))
* [`⊆`](@ref ⊆(::LazySet, ::LazySet))
* [`linear_combination`](@ref linear_combination(::LazySet, ::LazySet))
* [`minkowski_difference`](@ref minkowski_difference(::LazySet, ::LazySet))
* [`minkowski_sum`](@ref minkowski_sum(::LazySet, ::LazySet))

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`concretize`](@ref concretize(::LazySet))
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`eltype`](@ref eltype(::LazySet))
* [`extrema`](@ref extrema(::LazySet))
* [`extrema`](@ref extrema(::LazySet, ::Int))
* [`isoperation`](@ref isoperation(::LazySet))
* [`ispolyhedral`](@ref ispolyhedral(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`exact_sum`](@ref exact_sum(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))

Inherited from [`LazySet`](@ref) but does not apply:
* [`triangulate`](@ref triangulate(::LazySet))
