```@meta
CurrentModule = LazySets.UniverseModule
```

# [Universe](@id def_Universe)

```@docs
Universe
```

## Operations

```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
an_element(::LazySet)
```
```@meta
CurrentModule = LazySets.UniverseModule
```
```@docs
an_element(::Universe)
constrained_dimensions(::Universe)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
rand(::Type{LazySet})
```
```@meta
CurrentModule = LazySets.UniverseModule
```
```@docs
rand(::Type{Universe})
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
ρ(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets.UniverseModule
```
```@docs
ρ(::AbstractVector, ::Universe)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
σ(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets.UniverseModule
```
```@docs
σ(::AbstractVector, ::Universe)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:

* [`area`](@ref area(::LazySet))
```@meta
CurrentModule = LazySets
```
* [`chebyshev_center_radius`](@ref chebyshev_center_radius(::LazySet))
```@meta
CurrentModule = LazySets.API
```
* [`complement`](@ref complement(::LazySet))
* [`constraints_list`](@ref constraints_list(::LazySet))
* [`constraints`](@ref constraints(::LazySet))
* `copy(::Universe)`
* [`diameter`](@ref diameter(::LazySet, ::Real=Inf))
* [`dim`](@ref dim(::LazySet))
* [`isbounded`](@ref isbounded(::LazySet))
* [`isboundedtype`](@ref isboundedtype(::Type{LazySet}))
* [`isempty`](@ref isempty(::LazySet))
* [`isoperationtype`](@ref isoperationtype(::Type{LazySet}))
* [`isuniversal`](@ref isuniversal(::LazySet, ::Bool=false))
* [`norm`](@ref norm(::LazySet, ::Real=Inf))
* [`radius`](@ref radius(::LazySet, ::Real=Inf))
* [`rectify`](@ref rectify(::LazySet))
* [`reflect`](@ref reflect(::LazySet))
* [`distance`](@ref distance(::AbstractVector, ::LazySet))
* [`in`](@ref in(::AbstractVector, ::LazySet))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::LazySet))
* [`permute`](@ref permute(::LazySet, ::AbstractVector{Int}))
* [`project`](@ref project(::LazySet, ::AbstractVector{Int}))
* [`rationalize`](@ref rationalize(::LazySet))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`scale!`](@ref scale!(::Real, ::LazySet))
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
* [`translate!`](@ref translate!(::LazySet, ::AbstractVector))
* [`cartesian_product`](@ref cartesian_product(::LazySet, ::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet, ::LazySet))
* [`distance`](@ref distance(::LazySet, ::LazySet; ::Real=2.0))
* [`intersection`](@ref intersection(::LazySet, ::LazySet))

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`concretize`](@ref concretize(::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet))
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`eltype`](@ref eltype(::LazySet))
* [`isoperation`](@ref isoperation(::LazySet))
* [`ispolytopic`](@ref ispolytopic(::LazySet))
* [`ispolytopictype`](@ref ispolytopictype(::Type{LazySet}))
* [`polyhedron`](@ref polyhedron(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`tosimplehrep`](@ref tosimplehrep(::LazySet))
* [`triangulate`](@ref triangulate(::LazySet))
* [`triangulate_faces`](@ref triangulate_faces(::LazySet))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`sample`](@ref sample(::LazySet, ::Int=1))
* [`isapprox`](@ref isapprox(::LazySet, ::LazySet))
* [`isdisjoint`](@ref isdisjoint(::LazySet, ::LazySet))
* [`issubset`](@ref issubset(::LazySet, ::LazySet))
* [`minkowski_difference`](@ref minkowski_difference(::LazySet, ::LazySet))
* [`exact_sum`](@ref exact_sum(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))

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
* [`vertices`](@ref vertices(::AbstractPolyhedron))
* [`minkowski_sum`](@ref minkowski_sum(::AbstractPolyhedron, ::AbstractPolyhedron))
