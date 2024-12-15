```@meta
CurrentModule = LazySets.LineSegmentModule
```

# [Line segment (LineSegment)](@id def_LineSegment)

```@docs
LineSegment
an_element(::LineSegment)
center(::LineSegment)
constraints_list(::LineSegment)
dim(::LineSegment)
generators(::LineSegment)
genmat(::LineSegment)
halfspace_left(::LineSegment)
halfspace_right(::LineSegment)
ngens(::LineSegment)
rand(::Type{LineSegment})
vertices_list(::LineSegment)
∈(::AbstractVector, ::LineSegment)
σ(::AbstractVector, ::LineSegment)
ρ(::AbstractVector, ::LineSegment)
translate(::LineSegment, ::AbstractVector)
intersection(::LineSegment, ::LineSegment)
isdisjoint(::LineSegment, ::LineSegment)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:
* [`isoperationtype`](@ref isoperationtype(::Type{LazySet}))

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`area`](@ref area(::LazySet))
* [`complement`](@ref complement(::LazySet))
* [`concretize`](@ref concretize(::LazySet))
* [`constraints`](@ref constraints(::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`eltype`](@ref eltype(::LazySet))
* [`high`](@ref high(::LazySet))
* [`isoperation`](@ref isoperation(::LazySet))
* [`low`](@ref low(::LazySet))
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`reflect`](@ref reflect(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`surface`](@ref surface(::LazySet))
* [`vertices`](@ref vertices(::LazySet))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`sample`](@ref sample(::LazySet, ::Int))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`≈`](@ref ≈(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))

Inherited from [`ConvexSet`](@ref):
* [`high`](@ref high(::ConvexSet, ::Int))
* [`low`](@ref low(::ConvexSet, ::Int))
* [`linear_combination`](@ref linear_combination(::ConvexSet, ::ConvexSet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`ispolyhedral`](@ref ispolyhedral(::AbstractPolyhedron))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isboundedtype`](@ref isboundedtype(::Type{AbstractPolytope}))
* [`isconvextype`](@ref isconvextype(::Type{AbstractPolytope}))
* [`volume`](@ref volume(::AbstractPolytope))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`center`](@ref center(::AbstractCentrallySymmetricPolytope, ::Int))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetricPolytope))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetricPolytope, ::Int))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractCentrallySymmetricPolytope, ::Bool=false))

Inherited from [`AbstractZonotope`](@ref):
* [`order`](@ref order(::AbstractZonotope))
* [`reflect`](@ref reflect(::AbstractZonotope))
* [`togrep`](@ref togrep(::AbstractZonotope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractZonotope))
* [`cartesian_product`](@ref cartesian_product(::AbstractZonotope, ::AbstractZonotope))
* [`minkowski_difference`](@ref minkowski_difference(::AbstractZonotope, ::AbstractZonotope))
* [`minkowski_sum`](@ref minkowski_sum(::AbstractZonotope, ::AbstractZonotope))
