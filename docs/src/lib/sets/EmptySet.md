```@meta
CurrentModule = LazySets.EmptySetModule
```

# [Empty set (EmptySet)](@id def_EmptySet)

```@docs
EmptySet
∅
```

## Operations

```@docs
chebyshev_center_radius(::EmptySet; kwargs...)
complement(::EmptySet{N}) where {N}
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
* [`ρ`](@ref ρ(::AbstractVector, ::LazySet))
* [`σ`](@ref σ(::AbstractVector, ::LazySet))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`scale!`](@ref scale!(::Real, ::LazySet))
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
* [`translate!`](@ref translate!(::LazySet, ::AbstractVector))
* [`convex_hull`](@ref convex_hull(::LazySet, ::LazySet))
* [`intersection`](@ref isequivalent(::LazySet, ::LazySet))
* [`isdisjoint`](@ref isdisjoint(::LazySet, ::LazySet))
* [`⊆`](@ref ⊆(::LazySet, ::LazySet))

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
* [`cartesian_product`](@ref cartesian_product(::LazySet, ::LazySet))
* [`≈`](@ref ≈(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))
