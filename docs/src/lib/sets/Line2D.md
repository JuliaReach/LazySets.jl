```@meta
CurrentModule = LazySets
```

# [Line2D](@id def_Line2D)

```@docs
Line2D
```

## Operations

```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
an_element(::LazySet)
```
```@meta
CurrentModule = LazySets.Line2DModule
```
```@docs
an_element(::Line2D)
constrained_dimensions(::Line2D)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
isuniversal(::LazySet, ::Bool=false)
```
```@meta
CurrentModule = LazySets.Line2DModule
```
```@docs
isuniversal(::Line2D, ::Bool=false)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
rand(::Type{LazySet})
```
```@meta
CurrentModule = LazySets.Line2DModule
```
```@docs
rand(::Type{Line2D})
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
in(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets.Line2DModule
```
```@docs
in(::AbstractVector, ::Line2D)
project(::AbstractVector, ::Line2D)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
translate(::LazySet, ::AbstractVector)
```
```@meta
CurrentModule = LazySets.Line2DModule
```
```@docs
translate(::Line2D, ::AbstractVector)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
intersection(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets.Line2DModule
```
```@docs
intersection(::Line2D, ::Line2D)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:
* [`constraints_list`](@ref constraints_list(::LazySet))
* [`dim`](@ref dim(::LazySet))
* [`isbounded`](@ref isbounded(::LazySet))
* [`isempty`](@ref isempty(::LazySet))
* [`isoperationtype`](@ref isoperationtype(::Type{LazySet}))
* [`project`](@ref project(::LazySet, ::AbstractVector{Int}))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`scale!`](@ref scale!(::Real, ::LazySet))
* [`σ`](@ref σ(::AbstractVector, ::LazySet))
* [`isdisjoint`](@ref isdisjoint(::LazySet, ::LazySet))

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
* [`isboundedtype`](@ref isboundedtype(::Type{LazySet}))
* [`isconvex`](@ref isconvex(::LazySet))
* [`isoperation`](@ref isoperation(::LazySet))
* [`ispolytopic`](@ref ispolytopic(::LazySet))
* [`ispolytopictype`](@ref ispolytopictype(::Type{LazySet}))
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
* [`volume`](@ref volume(::LazySet))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::LazySet))
* [`sample`](@ref sample(::LazySet, ::Int=1))
* [`ρ`](@ref ρ(::AbstractVector, ::LazySet))
* [`translate!`](@ref translate!(::LazySet, ::AbstractVector))
* [`cartesian_product`](@ref cartesian_product(::LazySet, ::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet, ::LazySet))
* [`exact_sum`](@ref exact_sum(::LazySet, ::LazySet))
* [`isapprox`](@ref isapprox(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))
* [`minkowski_difference`](@ref minkowski_difference(::LazySet, ::LazySet))

Inherited from [`ConvexSet`](@ref):
* [`isconvextype`](@ref isconvextype(::Type{ConvexSet}))
* [`linear_combination`](@ref linear_combination(::ConvexSet, ::ConvexSet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`extrema`](@ref extrema(::AbstractPolyhedron))
* [`extrema`](@ref extrema(::AbstractPolyhedron, ::Int))
* [`high`](@ref high(::AbstractPolyhedron))
* [`high`](@ref high(::AbstractPolyhedron, ::Int))
* [`ispolyhedral`](@ref ispolyhedral(::AbstractPolyhedron))
* [`ispolyhedraltype`](@ref ispolyhedraltype(::Type{AbstractPolyhedron}))
* [`low`](@ref low(::AbstractPolyhedron))
* [`low`](@ref low(::AbstractPolyhedron, ::Int))
* [`vertices_list`](@ref vertices_list(::AbstractPolyhedron))
* [`issubset`](@ref issubset(::LazySet, ::AbstractPolyhedron))
* [`minkowski_sum`](@ref minkowski_sum(::AbstractPolyhedron, ::AbstractPolyhedron))
