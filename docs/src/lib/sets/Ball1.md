```@meta
CurrentModule = LazySets.Ball1Module
```

# [Manhattan-norm ball (Ball1)](@id def_Ball1)

```@docs
Ball1
```

## Operations

```@docs
center(::Ball1)
constraints_list(::Ball1)
rand(::Type{Ball1})
reflect(::Ball1)
vertices_list(::Ball1)
∈(::AbstractVector, ::Ball1, ::Bool=false)
ρ(::AbstractVector, ::Ball1)
σ(::AbstractVector, ::Ball1)
translate!(::Ball1, ::AbstractVector)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:
* [`high`](@ref high(::LazySet))
* [`high`](@ref high(::LazySet, ::Int))
* [`isoperationtype`](@ref isoperationtype(::Type{LazySet}))
* [`low`](@ref low(::LazySet))
* [`low`](@ref low(::LazySet, ::Int))
* [`project`](@ref project(::LazySet, ::AbstractVector{Int}))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`minkowski_sum`](@ref minkowski_sum(::LazySet, ::LazySet))
```@meta
CurrentModule = LazySets
```
* [`ball_norm`](@ref ball_norm(::LazySet))
* [`radius_ball`](@ref radius_ball(::LazySet))

Inherited from [`LazySet`](@ref):
* [`area`](@ref area(::LazySet))
* [`complement`](@ref complement(::LazySet))
* [`concretize`](@ref concretize(::LazySet))
* [`constraints`](@ref constraints(::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`eltype`](@ref eltype(::LazySet))
* [`rectify`](@ref rectify(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`surface`](@ref surface(::LazySet))
* [`vertices`](@ref vertices(::LazySet))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::LazySet)
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`sample`](@ref sample(::LazySet, ::Int))
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
* [`cartesian_product`](@ref cartesian_product(::LazySet, ::LazySet))
* [`≈`](@ref ≈(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))
* [`minkowski_difference`](@ref minkowski_difference(::LazySet, ::LazySet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`is_polyhedral`](@ref is_polyhedral(::AbstractPolyhedron))
* [`intersection`](@ref intersection(::AbstractPolyhedron, ::AbstractPolyhedron))
* [`isdisjoint`](@ref isdisjoint(::AbstractPolyhedron, ::AbstractPolyhedron))
* [`⊆`](@ref ⊆(::AbstractPolyhedron, ::AbstractPolyhedron))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
```@meta
CurrentModule = LazySets.API
```
* [`isboundedtype`](@ref isboundedtype(::Type{LazySet}))
* [`isconvextype`](@ref isconvextype(::Type{LazySet}))
```@meta
CurrentModule = LazySets
```
* [`isoperation`](@ref isoperation(::AbstractPolytope))
* [`volume`](@ref volume(::AbstractPolytope))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope))
* [`center`](@ref center(::AbstractCentrallySymmetricPolytope, ::Int))
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetricPolytope))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetricPolytope, ::Int))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractCentrallySymmetricPolytope{N}, ::Bool=false) where {N})
