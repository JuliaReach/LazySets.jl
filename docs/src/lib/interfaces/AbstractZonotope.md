```@contents
Pages = ["AbstractZonotope.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# [Zonotopes (AbstractZonotope)](@id def_AbstractZonotope)

A zonotope is a specific centrally symmetric polytope characterized by a
center and a collection of generators.

```@docs
AbstractZonotope
```

This interface requires to implement the following functions:

```@docs
generators(::AbstractZonotope)
genmat(::AbstractZonotope)
```

This interface defines the following functions:

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
constraints_list(::AbstractZonotope)
constraints_list(::AbstractZonotope{<:AbstractFloat}; ::Bool=true)
ngens(::AbstractZonotope)
order(::AbstractZonotope)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
reflect(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
reflect(::AbstractZonotope)
remove_redundant_generators(::AbstractZonotope)
togrep(::AbstractZonotope)
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
vertices_list(::AbstractZonotope)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
in(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
in(::AbstractVector, ::AbstractZonotope)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
linear_map(::AbstractMatrix, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
linear_map(::AbstractMatrix, ::AbstractZonotope)
reduce_order(::AbstractZonotope, ::Real, ::AbstractReductionMethod=GIR05())
split(::AbstractZonotope, ::Int)
split(::AbstractZonotope, ::AbstractVector{Int}, ::AbstractVector{Int})
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
ρ(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
ρ(::AbstractVector, ::AbstractZonotope)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
σ(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
σ(::AbstractVector, ::AbstractZonotope)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
isdisjoint(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
isdisjoint(::AbstractZonotope, ::AbstractZonotope, ::Bool=false)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
minkowski_difference(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
minkowski_difference(::AbstractZonotope, ::AbstractZonotope)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
minkowski_sum(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
minkowski_sum(::AbstractZonotope, ::AbstractZonotope)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
norm(::LazySet, ::Real=Inf)
```
```@meta
CurrentModule = LazySets
```
```@docs
norm(::AbstractZonotope, ::Real=Inf)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:

* [`project`](@ref project(::LazySet, ::AbstractVector{Int}))
* [`cartesian_product`](@ref cartesian_product(::LazySet, ::LazySet))

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
* [`ispolytopic`](@ref ispolytopic(::LazySet))
* [`polyhedron`](@ref polyhedron(::LazySet))
* [`radius`](@ref radius(::LazySet, ::Real=Inf))
* [`rationalize`](@ref rationalize(::LazySet))
* [`rectify`](@ref rectify(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`tosimplehrep`](@ref tosimplehrep(::LazySet))
* [`triangulate`](@ref triangulate(::LazySet))
* [`triangulate_faces`](@ref triangulate_faces(::LazySet))
* [`vertices`](@ref vertices(::LazySet))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`sample`](@ref sample(::LazySet, ::Int=1))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
* [`convex_hull`](@ref convex_hull(::LazySet, ::LazySet))
* [`exact_sum`](@ref exact_sum(::LazySet, ::LazySet))
* [`isapprox`](@ref isapprox(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))

Inherited from [`ConvexSet`](@ref):
* [`isconvextype`](@ref isconvextype(::Type{LazySet}))
* [`linear_combination`](@ref linear_combination(::ConvexSet, ::ConvexSet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`high`](@ref high(::AbstractPolyhedron))
* [`high`](@ref high(::AbstractPolyhedron, ::Int))
* [`ispolyhedral`](@ref ispolyhedral(::AbstractPolyhedron))
* [`ispolyhedraltype`](@ref ispolyhedraltype(::Type{AbstractPolyhedron}))
* [`low`](@ref low(::AbstractPolyhedron))
* [`low`](@ref low(::AbstractPolyhedron, ::Int))
* [`intersection`](@ref intersection(::AbstractPolyhedron{N}, ::AbstractPolyhedron{N}) where {N})

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isboundedtype`](@ref isboundedtype(::Type{AbstractPolytope}))
* [`volume`](@ref volume(::AbstractPolytope))
* [`⊆`](@ref ⊆(::AbstractPolytope, ::LazySet, ::Bool=false))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope))
* [`center`](@ref center(::AbstractCentrallySymmetricPolytope, ::Int))
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetricPolytope))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetricPolytope, ::Int))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractCentrallySymmetricPolytope, ::Bool=false))

## Internal methods

```@docs
generators_fallback(::AbstractZonotope)
genmat_fallback(::AbstractZonotope)
_l1_norm(::AbstractZonotope)
```

## Order-reduction methods

```@docs
LazySets.AbstractReductionMethod
LazySets.ASB10
LazySets.COMB03
LazySets.GIR05
LazySets.SRMB16
```

## Implementations

* [Zonotope](@ref def_Zonotope)
* [Line segment (LineSegment)](@ref def_LineSegment)
