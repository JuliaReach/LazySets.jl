```@contents
Pages = ["AbstractAffineMap.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# [Affine maps (AbstractAffineMap)](@id def_AbstractAffineMap)

An affine map consists of a linear map and a translation.

```@docs
AbstractAffineMap
```

This interface requires to implement the following functions:

```@docs
matrix(::AbstractAffineMap)
vector(::AbstractAffineMap)
set(::AbstractAffineMap)
```

This interface defines the following functions:

```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
an_element(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
an_element(::AbstractAffineMap)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
center(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
center(::AbstractAffineMap)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
constraints_list(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
constraints_list(::AbstractAffineMap)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
isbounded(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
isbounded(::AbstractAffineMap)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
isempty(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
isempty(::AbstractAffineMap)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
∈(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
∈(::AbstractVector, ::AbstractAffineMap)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
vertices_list(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
vertices_list(::AbstractAffineMap)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:

* [`dim`](@ref dim(::LazySet))
* [`isconvextype`](@ref isconvextype(::Type{<:LazySet}))
* [`isoperationtype`](@ref isoperationtype(::Type{<:LazySet}))
* [`ispolyhedral`](@ref ispolyhedral(::LazySet))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::LazySet))
* [`ρ`](@ref ρ(::AbstractVector, ::LazySet))
* [`σ`](@ref σ(::AbstractVector, ::LazySet))

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`area`](@ref area(::LazySet))
* [`chebyshev_center_radius`](@ref chebyshev_center_radius(::LazySet))
* [`complement`](@ref complement(::LazySet))
* [`concretize`](@ref concretize(::LazySet))
* [`constraints`](@ref constraints(::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet))
* `copy(::Type{LazySet})`
* [`diameter`](@ref diameter(::LazySet, ::Real=Inf))
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`eltype`](@ref eltype(::LazySet))
* [`extrema`](@ref extrema(::LazySet))
* [`extrema`](@ref extrema(::LazySet, ::Int))
* [`high`](@ref high(::LazySet))
* [`high`](@ref high(::LazySet, ::Int))
* [`isboundedtype`](@ref isboundedtype(::Type{LazySet}))
* [`isoperation`](@ref isoperation(::LazySet))
* [`low`](@ref low(::LazySet))
* [`low`](@ref low(::LazySet, ::Int))
* [`norm`](@ref norm(::LazySet, ::Real=Inf))
* [`polyhedron`](@ref polyhedron(::LazySet))
* [`radius`](@ref radius(::LazySet, ::Real=Inf))
* [`rationalize`](@ref rationalize(::LazySet))
* [`rectify`](@ref rectify(::LazySet))
* [`reflect`](@ref reflect(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`tosimplehrep`](@ref tosimplehrep(::LazySet))
* [`triangulate`](@ref triangulate(::LazySet))
* [`triangulate_faces`](@ref triangulate_faces(::LazySet))
* [`vertices`](@ref vertices(::LazySet))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`project`](@ref project(::LazySet, ::AbstractVector{Int}))
* [`sample`](@ref sample(::LazySet, ::Int=1))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
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

## Implementations

* [Affine map (AffineMap)](@ref def_AffineMap)
* [Exponential map (ExponentialMap)](@ref def_ExponentialMap)
* [Linear map (LinearMap)](@ref def_LinearMap)
* [Reset map (ResetMap)](@ref def_ResetMap)
* [Translation](@ref def_Translation)
