```@contents
Pages = ["AbstractCentrallySymmetricPolytope.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# [Centrally symmetric polytopes (AbstractCentrallySymmetricPolytope)](@id def_AbstractCentrallySymmetricPolytope)

A centrally symmetric polytope is a combination of two other interfaces:
[Centrally symmetric sets](@ref def_AbstractCentrallySymmetric) and
[Polytope](@ref def_AbstractPolytope).

```@docs
AbstractCentrallySymmetricPolytope
```

This interface requires to implement the required functions of both the
`AbstractCentrallySymmetric` and `AbstractPolytope` interfaces:

```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
center(::LazySet)
constraints_list(::LazySet)
vertices_list(::LazySet)
```

```@meta
CurrentModule = LazySets
```

This interface shares the following functions with
[`AbstractCentrallySymmetric`](@ref):

```@docs
an_element(::AbstractCentrallySymmetricPolytope)
extrema(::AbstractCentrallySymmetricPolytope)
extrema(::AbstractCentrallySymmetricPolytope, ::Int)
isuniversal(::AbstractCentrallySymmetricPolytope, ::Bool=false)
```

Undocumented implementations shared with [`AbstractCentrallySymmetric`](@ref):

```@meta
CurrentModule = LazySets.API
```

* [`center`](@ref center(::LazySet, ::Int))
* [`dim`](@ref dim(::LazySet))
* [`isempty`](@ref isempty(::LazySet))

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
* [`isoperation`](@ref isoperation(::LazySet))
* [`norm`](@ref norm(::LazySet, ::Real=Inf))
* [`radius`](@ref radius(::LazySet, ::Real=Inf))
* [`rationalize`](@ref rationalize(::LazySet))
* [`rectify`](@ref rectify(::LazySet))
* [`reflect`](@ref reflect(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`surface`](@ref surface(::LazySet))
* [`triangulate`](@ref triangulate(::LazySet))
* [`vertices`](@ref vertices(::LazySet))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::LazySet))
* [`sample`](@ref sample(::LazySet, ::Int=1))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`ρ`](@ref ρ(::AbstractVector, ::LazySet))
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
* [`cartesian_product`](@ref cartesian_product(::LazySet, ::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet, ::LazySet))
* [`exact_sum`](@ref exact_sum(::LazySet, ::LazySet))
* [`≈`](@ref ≈(::LazySet, ::LazySet))
* [`isdisjoint`](@ref isdisjoint(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))
* [`minkowski_difference`](@ref minkowski_difference(::LazySet, ::LazySet))

Inherited from [`ConvexSet`](@ref):
* [`isconvextype`](@ref isconvextype(::Type{LazySet}))
* [`linear_combination`](@ref linear_combination(::ConvexSet, ::ConvexSet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`high`](@ref high(::AbstractPolyhedron))
* [`high`](@ref high(::AbstractPolyhedron, ::Int))
* [`ispolyhedral`](@ref ispolyhedral(::AbstractPolyhedron))
* [`low`](@ref low(::AbstractPolyhedron))
* [`low`](@ref low(::AbstractPolyhedron, ::Int))
* [`∈`](@ref ∈(::AbstractVector, ::AbstractPolyhedron))
* [`project`](@ref project(::AbstractPolyhedron, ::AbstractVector{Int}))
* [`minkowski_sum`](@ref minkowski_sum(::AbstractPolyhedron, ::AbstractPolyhedron))
* [`intersection`](@ref intersection(::AbstractPolyhedron{N}, ::AbstractPolyhedron{N}) where {N})

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isboundedtype`](@ref isboundedtype(::Type{AbstractPolytope}))
* [`volume`](@ref volume(::AbstractPolytope))
* [`⊆`](@ref ⊆(::AbstractPolytope, ::LazySet, ::Bool=false))

## Implementations

* [Manhattan-norm ball (Ball1)](@ref def_Ball1)
