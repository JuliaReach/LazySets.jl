```@contents
Pages = ["AbstractPolytope.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# [Polytopes (AbstractPolytope)](@id def_AbstractPolytope)

A polytope is a bounded set with finitely many vertices (*V-representation*)
resp. facets (*H-representation*).
Note that there is a special interface combination
[Centrally symmetric polytope](@ref def_AbstractCentrallySymmetricPolytope).

```@docs
AbstractPolytope
```

This interface requires to implement the following function:

```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
vertices_list(::LazySet)
```

This interface defines the following functions:

```@docs; canonical=false
isempty(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
isempty(::AbstractPolytope)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
isuniversal(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
isuniversal(::AbstractPolytope, ::Bool=false)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
volume(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
volume(::AbstractPolytope)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
⊆(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
⊆(::AbstractPolytope, ::LazySet, ::Bool=false)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:

* [`isbounded`](@ref isbounded(::LazySet))
* [`isboundedtype`](@ref isboundedtype(::Type{LazySet}))

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
* [`⊆`](@ref ⊆(::LazySet, ::LazySet))
* [`minkowski_difference`](@ref minkowski_difference(::LazySet, ::LazySet))
* [`minkowski_sum`](@ref minkowski_sum(::LazySet, ::LazySet))

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

Some common functions implemented by several subtypes:

```@docs
remove_redundant_vertices(::AbstractPolytope)
remove_redundant_vertices!(::AbstractPolytope)
```

## Implementations

* [Polytope in constraint representation (HPolytope)](@ref def_HPolytope)
* [Polytope in vertex representation (VPolytope)](@ref def_VPolytope)
