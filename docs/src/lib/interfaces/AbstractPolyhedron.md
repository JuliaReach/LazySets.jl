```@contents
Pages = ["AbstractPolyhedron.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# [Polyhedra (AbstractPolyhedron)](@id def_AbstractPolyhedron)

A polyhedron has finitely many facets (*H-representation*) and is not
necessarily bounded.

```@docs
AbstractPolyhedron
```

This interface requires to implement the following function:

```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
constraints_list(::LazySet)
```
```@meta
CurrentModule = LazySets
```

This interface defines the following functions:

```@docs
an_element(::AbstractPolyhedron)
constrained_dimensions(::AbstractPolyhedron)
isbounded(::AbstractPolyhedron)
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
isuniversal(::AbstractPolyhedron, ::Bool=false)
vertices_list(::AbstractPolyhedron)
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
in(::AbstractVector, ::AbstractPolyhedron)
project(::AbstractPolyhedron, ::AbstractVector{Int})
intersection(::AbstractPolyhedron{N}, ::AbstractPolyhedron{N}) where {N}
minkowski_sum(::AbstractPolyhedron, ::AbstractPolyhedron)
LazySets._isbounded_stiemke
LazySets._linear_map_polyhedron
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:

* [`extrema`](@ref extrema(::LazySet))
* [`extrema`](@ref extrema(::LazySet, ::Int))
* [`high`](@ref high(::LazySet))
* [`high`](@ref high(::LazySet, ::Int))
* [`ispolyhedral`](@ref ispolyhedral(::LazySet))
* [`ispolyhedraltype`](@ref ispolyhedraltype(::Type{LazySet}))
* [`low`](@ref low(::LazySet))
* [`low`](@ref low(::LazySet, ::Int))
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
* [`isempty`](@ref isempty(::LazySet))
* [`isoperation`](@ref isoperation(::LazySet))
* [`ispolytopic`](@ref ispolytopic(::LazySet))
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
* [`isapprox`](@ref isapprox(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))
* [`issubset`](@ref issubset(::LazySet, ::LazySet))
* [`minkowski_difference`](@ref minkowski_difference(::LazySet, ::LazySet))

Inherited from [`ConvexSet`](@ref):
* [`isconvextype`](@ref isconvextype(::Type{LazySet}))
* [`linear_combination`](@ref linear_combination(::ConvexSet, ::ConvexSet))

Some common functions implemented by several subtypes:

```@docs
addconstraint!(::AbstractPolyhedron, ::HalfSpace)
ishyperplanar(::AbstractPolyhedron)
```

Plotting polyhedra is available too:

```@docs
plot_recipe(::AbstractPolyhedron{N}, ::Any=zero(N)) where {N}
```

## Implementations

* [Half-space (HalfSpace)](@ref def_HalfSpace)
* [Polyhedron in constraint representation (HPolyhedron)](@ref def_HPolyhedron)
* [Hyperplane](@ref def_Hyperplane)
* [Line2D](@ref def_Line2D)
* [Line](@ref def_Line)
* [Universe](@ref def_Universe)
* [Star](@ref def_Star)
