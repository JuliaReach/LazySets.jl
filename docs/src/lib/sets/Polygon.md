```@meta
CurrentModule = LazySets.PolygonModule
```

# [Non-convex polygon (Polygon)](@id def_Polygon)

```@docs
Polygon
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:

* [`convex_hull`](@ref convex_hull(::LazySet))
* [`dim`](@ref dim(::LazySet))
* [`isbounded`](@ref isbounded(::LazySet))
* [`isboundedtype`](@ref isboundedtype(::Type{LazySet}))
* [`isconvextype`](@ref isconvextype(::Type{LazySet}))
* [`isempty`](@ref isempty(::LazySet))
* [`isoperationtype`](@ref isoperationtype(::Type{LazySet}))
* [`∈`](@ref ∈(::AbstractVector, ::LazySet))
* [`ρ`](@ref ρ(::AbstractVector, ::LazySet))
* [`σ`](@ref σ(::AbstractVector, ::LazySet))

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`an_element`](@ref an_element(::LazySet))
* [`area`](@ref area(::LazySet))
* [`center`](@ref center(::LazySet))
* [`complement`](@ref complement(::LazySet))
* [`concretize`](@ref concretize(::LazySet))
* [`constraints`](@ref constraints(::LazySet))
* `copy(::Type{LazySet})`
* [`diameter`](@ref diameter(::LazySet, ::Real=Inf))
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`eltype`](@ref eltype(::LazySet))
* [`extrema`](@ref extrema(::LazySet))
* [`extrema`](@ref extrema(::LazySet, ::Int))
* [`high`](@ref high(::LazySet))
* [`high`](@ref high(::LazySet, ::Int))
* [`isoperation`](@ref isoperation(::LazySet))
* [`ispolyhedral`](@ref ispolyhedral(::LazySet))
* [`low`](@ref low(::LazySet))
* [`low`](@ref low(::LazySet, ::Int))
* [`norm`](@ref norm(::LazySet, ::Real=Inf))
* [`radius`](@ref radius(::LazySet, ::Real=Inf))
* [`rectify`](@ref rectify(::LazySet))
* [`reflect`](@ref reflect(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`surface`](@ref surface(::LazySet))
* [`vertices`](@ref vertices(::LazySet))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::LazySet))
* [`project`](@ref project(::LazySet, ::AbstractVector{Int}))
* [`sample`](@ref sample(::LazySet, ::Int=1))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`scale!`](@ref scale!(::Real, ::LazySet))
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
* [`translate!`](@ref translate!(::LazySet, ::AbstractVector))
* [`cartesian_product`](@ref cartesian_product(::LazySet, ::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet, ::LazySet))
* [`exact_sum`](@ref exact_sum(::LazySet, ::LazySet))
* [`≈`](@ref ≈(::LazySet, ::LazySet))
* [`isdisjoint`](@ref isdisjoint(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))
* [`⊆`](@ref ⊆(::LazySet, ::LazySet))
* [`linear_combination`](@ref linear_combination(::LazySet, ::LazySet))
* [`minkowski_difference`](@ref minkowski_difference(::LazySet, ::LazySet))
* [`minkowski_sum`](@ref minkowski_sum(::LazySet, ::LazySet))
