```@meta
CurrentModule = LazySets.LineSegmentModule
```

# [Line segment (LineSegment)](@id def_LineSegment)

```@docs
LineSegment
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
an_element(::LazySet)
```
```@meta
CurrentModule = LazySets.LineSegmentModule
```
```@docs
an_element(::LineSegment)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
constraints_list(::LazySet)
```
```@meta
CurrentModule = LazySets.LineSegmentModule
```
```@docs
constraints_list(::LineSegment)
generators(::LineSegment)
genmat(::LineSegment)
halfspace_left(::LineSegment)
halfspace_right(::LineSegment)
ngens(::LineSegment)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
rand(::Type{LazySet})
```
```@meta
CurrentModule = LazySets.LineSegmentModule
```
```@docs
rand(::Type{LineSegment})
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
∈(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets.LineSegmentModule
```
```@docs
∈(::AbstractVector, ::LineSegment)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
σ(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets.LineSegmentModule
```
```@docs
σ(::AbstractVector, ::LineSegment)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
translate(::LazySet, ::AbstractVector)
```
```@meta
CurrentModule = LazySets.LineSegmentModule
```
```@docs
translate(::LineSegment, ::AbstractVector)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
intersection(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets.LineSegmentModule
```
```@docs
intersection(::LineSegment, ::LineSegment)
isdisjoint(::LineSegment, ::LineSegment, ::Bool=false)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:
* [`center`](@ref center(::LazySet))
* [`dim`](@ref dim(::LazySet))
* [`isoperationtype`](@ref isoperationtype(::Type{LazySet}))
* [`vertices_list`](@ref vertices_list(::LazySet))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`scale!`](@ref scale!(::Real, ::LazySet))
* [`ρ`](@ref ρ(::AbstractVector, ::LazySet))

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
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`eltype`](@ref eltype(::LazySet))
* [`isoperation`](@ref isoperation(::LazySet))
* [`norm`](@ref norm(::LazySet, ::Real))
* [`polyhedron`](@ref polyhedron(::LazySet))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`rationalize`](@ref rationalize(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`tosimplehrep`](@ref tosimplehrep(::LazySet))
* [`triangulate`](@ref triangulate(::LazySet))
* [`triangulate_faces`](@ref triangulate_faces(::LazySet))
* [`vertices`](@ref vertices(::LazySet))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`sample`](@ref sample(::LazySet, ::Int))
* [`convex_hull`](@ref convex_hull(::LazySet, ::LazySet))
* [`≈`](@ref ≈(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))

Inherited from [`ConvexSet`](@ref):
* [`isconvextype`](@ref isconvextype(::Type{ConvexSet}))
* [`linear_combination`](@ref linear_combination(::ConvexSet, ::ConvexSet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`high`](@ref high(::AbstractPolyhedron))
* [`high`](@ref high(::AbstractPolyhedron, ::Int))
* [`ispolyhedral`](@ref ispolyhedral(::AbstractPolyhedron))
* [`low`](@ref low(::AbstractPolyhedron))
* [`low`](@ref low(::AbstractPolyhedron, ::Int))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isboundedtype`](@ref isboundedtype(::Type{LazySet}))
* [`volume`](@ref volume(::AbstractPolytope))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`center`](@ref center(::AbstractCentrallySymmetric, ::Int))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetricPolytope))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetricPolytope, ::Int))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractCentrallySymmetricPolytope, ::Bool=false))

Inherited from [`AbstractZonotope`](@ref):
* [`order`](@ref order(::AbstractZonotope))
* [`reflect`](@ref reflect(::AbstractZonotope))
* [`remove_redundant_generators`](@ref remove_redundant_generators(::AbstractZonotope))
* [`togrep`](@ref togrep(::AbstractZonotope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractZonotope))
* [`project`](@ref project(::AbstractZonotope, ::AbstractVector{Int}))
* [`reduce_order`](@ref reduce_order(::AbstractZonotope, ::Real, ::AbstractReductionMethod=GIR05()))
* [`split`](@ref split(::AbstractZonotope, ::Int))
* [`split`](@ref split(::AbstractZonotope, ::AbstractVector{Int}, ::AbstractVector{Int}))
* [`cartesian_product`](@ref cartesian_product(::AbstractZonotope, ::AbstractZonotope))
* [`minkowski_difference`](@ref minkowski_difference(::AbstractZonotope, ::AbstractZonotope))
* [`minkowski_sum`](@ref minkowski_sum(::AbstractZonotope, ::AbstractZonotope))
