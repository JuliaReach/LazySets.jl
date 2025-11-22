```@meta
CurrentModule = LazySets.HPolyhedronModule
```

# [Polyhedron in constraint representation (HPolyhedron)](@id def_HPolyhedron)

```@docs
HPolyhedron
```

## Conversion

```julia
convert(::Type{HPolyhedron}, ::LazySet)
convert(::Type{HPolyhedron}, ::Polyhedra.HRep)
```

## Operations

```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
rand(::Type{LazySet})
```
```@meta
CurrentModule = LazySets.HPolyhedronModule
```
```@docs
rand(::Type{HPolyhedron})
```

The following methods are shared between `HPolytope` and `HPolyhedron`.

```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
dim(::LazySet)
```
```@meta
CurrentModule = LazySets.HPolyhedronModule
```
```@docs
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
```@docs
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
```@docs
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
```@docs
translate(::HPoly, ::AbstractVector)
convex_hull(::HPoly, ::HPoly)
```

Undocumented implementations:
* [`constraints_list`](@ref constraints_list(::LazySet))
* [`isempty`](@ref isempty(::LazySet))
`ishyperplanar(::HPolyhedron)`
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
* [`diameter`](@ref diameter(::LazySet, ::Real=Inf))
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`eltype`](@ref eltype(::LazySet))
* [`isboundedtype`](@ref isboundedtype(::Type{LazySet}))
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
* [`volume`](@ref volume(::LazySet))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::LazySet))
* [`sample`](@ref sample(::LazySet, ::Int=1))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`cartesian_product`](@ref cartesian_product(::LazySet, ::LazySet))
* [`difference`](@ref difference(::LazySet, ::LazySet))
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
* [`an_element`](@ref an_element(::AbstractPolyhedron))
* [`constrained_dimensions`](@ref constrained_dimensions(::AbstractPolyhedron))
* [`extrema`](@ref extrema(::AbstractPolyhedron))
* [`extrema`](@ref extrema(::AbstractPolyhedron, ::Int))
* [`high`](@ref high(::AbstractPolyhedron))
* [`high`](@ref high(::AbstractPolyhedron, ::Int))
* [`isbounded`](@ref isbounded(::AbstractPolyhedron))
* [`ispolyhedral`](@ref ispolyhedral(::AbstractPolyhedron))
* [`ispolyhedraltype`](@ref ispolyhedraltype(::Type{AbstractPolyhedron}))
* [`isuniversal`](@ref isuniversal(::AbstractPolyhedron, ::Bool=false))
* [`low`](@ref low(::AbstractPolyhedron))
* [`low`](@ref low(::AbstractPolyhedron, ::Int))
* [`vertices_list`](@ref vertices_list(::AbstractPolyhedron, ::Bool=false))
* [`in`](@ref in(::AbstractVector, ::AbstractPolyhedron))
* [`project`](@ref project(::AbstractPolyhedron, ::AbstractVector{Int}))
* [`intersection`](@ref intersection(::AbstractPolyhedron, ::AbstractPolyhedron))
* [`isdisjoint`](@ref isdisjoint(::AbstractPolyhedron, ::AbstractPolyhedron))
* [`issubset`](@ref issubset(::LazySet, ::AbstractPolyhedron))
* [`minkowski_sum`](@ref minkowski_sum(::AbstractPolyhedron, ::AbstractPolyhedron))
