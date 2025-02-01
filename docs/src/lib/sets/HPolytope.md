```@meta
CurrentModule = LazySets
```

# [Polytope in constraint representation (HPolytope)](@id def_HPolytope)

[Convex polytopes](https://en.wikipedia.org/wiki/Polytope) are bounded polyhedra.
The type `HPolytope` represents polytopes.
While identical to [`HPolyhedron`](@ref) in its representation, `HPolytope`
instances are assumed to be bounded.

```@meta
CurrentModule = LazySets.HPolytopeModule
```

```@docs
HPolytope
```

## Conversion

```julia
convert(::Type{HPolytope}, ::LazySet)
convert(::Type{HPolytope}, ::Polyhedra.HRep)
```

## Operations

```@meta
CurrentModule = LazySets.HPolytopeModule
```

```@docs
isbounded(::HPolytope, ::Bool=true)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
rand(::Type{LazySet})
```
```@meta
CurrentModule = LazySets.HPolytopeModule
```
```@docs
rand(::Type{HPolytope})
vertices_list(::HPolytope)
```

```@meta
CurrentModule = LazySets.HPolyhedronModule
```

The following functionality is shared with [`HPolyhedron`](@ref).

```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
dim(::LazySet)
```
```@meta
CurrentModule = LazySets.HPolyhedronModule
```
```@docs; canonical=false
dim(::HPoly)
normalize(::HPoly{N}, p::Real=N(2)) where {N}
remove_redundant_constraints(::HPoly)
remove_redundant_constraints!(::HPoly)
tohrep(::HPoly)
tovrep(::HPoly)
addconstraint!(::HPoly, ::HalfSpace)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
ρ(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets.HPolyhedronModule
```
```@docs; canonical=false
ρ(::AbstractVector{M}, ::HPoly{N}) where {M, N}
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
σ(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets.HPolyhedronModule
```
```@docs; canonical=false
σ(::AbstractVector{M}, ::HPoly{N}) where {M, N}
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
translate(::LazySet, ::AbstractVector)
```
```@meta
CurrentModule = LazySets.HPolyhedronModule
```
```@docs; canonical=false
translate(::HPoly, ::AbstractVector)
convex_hull(::HPoly, ::HPoly)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:
* [`constraints_list`](@ref constraints_list(::LazySet))
* [`isempty`](@ref isempty(::LazySet))
* [`isoperationtype`](@ref isoperationtype(::Type{LazySet}))
* [`permute`](@ref permute(::LazySet, ::AbstractVector{Int}))
* [`minkowski_sum`](@ref minkowski_sum(::LazySet, ::LazySet))

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
* [`cartesian_product`](@ref cartesian_product(::LazySet, ::LazySet))
* [`exact_sum`](@ref exact_sum(::LazySet, ::LazySet))
* [`≈`](@ref ≈(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))
* [`minkowski_difference`](@ref minkowski_difference(::LazySet, ::LazySet))

Inherited from [`ConvexSet`](@ref):
* [`isconvextype`](@ref isconvextype(::Type{ConvexSet}))
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
* [`intersection`](@ref intersection(::AbstractPolyhedron, ::AbstractPolyhedron))
* [`isdisjoint`](@ref isdisjoint(::AbstractPolyhedron, ::AbstractPolyhedron))

Inherited from [`AbstractPolytope`](@ref):
* [`isboundedtype`](@ref isboundedtype(::Type{AbstractPolytope}))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope, ::Bool=false))
* [`volume`](@ref volume(::AbstractPolytope))
* [`⊆`](@ref ⊆(::AbstractPolytope, ::AbstractPolyhedron))
