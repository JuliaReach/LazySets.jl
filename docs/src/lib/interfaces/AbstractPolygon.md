```@contents
Pages = ["AbstractPolygon.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# [Polygons (AbstractPolygon)](@id def_AbstractPolygon)

A polygon is a two-dimensional polytope.

```@docs
AbstractPolygon
```

This interface defines the following functions (undocumented):

```@meta
CurrentModule = LazySets.API
```

* [`dim`](@ref dim(::LazySet))
* [`volume`](@ref volume(::LazySet))
* [`intersection`](@ref intersection(::LazySet, ::LazySet))

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`area`](@ref area(::LazySet))
* [`center`](@ref center(::LazySet))
* [`chebyshev_center_radius`](@ref chebyshev_center_radius(::LazySet))
* [`complement`](@ref complement(::LazySet))
* [`concretize`](@ref concretize(::LazySet))
* [`constraints`](@ref constraints(::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet))
* `copy(::Type{LazySet})`
* [`delaunay`](@ref delaunay(::LazySet))
* [`diameter`](@ref diameter(::LazySet, ::Real=Inf))
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`eltype`](@ref eltype(::LazySet))
* [`isoperation`](@ref isoperation(::LazySet))
* [`norm`](@ref norm(::LazySet, ::Real=Inf))
* [`polyhedron`](@ref polyhedron(::LazySet))
* [`radius`](@ref radius(::LazySet, ::Real=Inf))
* [`rationalize`](@ref rationalize(::LazySet))
* [`rectify`](@ref rectify(::LazySet))
* [`reflect`](@ref reflect(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`surface`](@ref surface(::LazySet))
* [`tosimplehrep`](@ref tosimplehrep(::LazySet))
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
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))
* [`minkowski_difference`](@ref minkowski_difference(::LazySet, ::LazySet))

Inherited from [`ConvexSet`](@ref):
* [`isconvextype`](@ref isconvextype(::Type{LazySet}))
* [`linear_combination`](@ref linear_combination(::ConvexSet, ::ConvexSet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`an_element`](@ref an_element(::AbstractPolyhedron))
* [`extrema`](@ref extrema(::AbstractPolyhedron))
* [`extrema`](@ref extrema(::AbstractPolyhedron, ::Int))
* [`high`](@ref high(::AbstractPolyhedron))
* [`high`](@ref high(::AbstractPolyhedron, ::Int))
* [`ispolyhedral`](@ref ispolyhedral(::AbstractPolyhedron))
* [`low`](@ref low(::AbstractPolyhedron))
* [`low`](@ref low(::AbstractPolyhedron, ::Int))
* [`∈`](@ref ∈(::AbstractVector, ::AbstractPolyhedron))
* [`project`](@ref project(::AbstractPolyhedron, ::AbstractVector{Int}))
* [`isdisjoint`](@ref isdisjoint(::AbstractPolyhedron, ::AbstractPolyhedron))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::LazySet))
* [`isboundedtype`](@ref isboundedtype(::Type{LazySet}))
* [`isempty`](@ref isempty(::LazySet, ::Bool=false))
* [`isuniversal`](@ref isuniversal(::LazySet, ::Bool=false))
* [`⊆`](@ref ⊆(::AbstractPolytope, ::LazySet, ::Bool=false))
* [`minkowski_sum`](@ref minkowski_sum(::AbstractPolytope, ::AbstractPolytope))

The following helper functions are used for sorting directions:

```@docs
LazySets.jump2pi
⪯(::AbstractVector, ::AbstractVector)
LazySets._leq_trig(::AbstractVector{N}, ::AbstractVector{N}) where {N<:AbstractFloat}
LazySets.quadrant(::AbstractVector)
```

## Implementations

* [Polygon in vertex representation (VPolygon)](@ref def_VPolygon)
