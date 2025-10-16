```@meta
CurrentModule = LazySets.VPolytopeModule
```

# [Polytope in vertex representation (VPolytope)](@id def_VPolytope)

```@docs
VPolytope
```

## Conversion

```@docs
convert(::Type{VPolytope}, ::LazySet)
```
```julia
convert(::Type{VPolytope}, ::Polyhedra.VRep)
```

## Operations

```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
constraints_list(::LazySet)
```
```@meta
CurrentModule = LazySets.VPolytopeModule
```
```@docs
constraints_list(::VPolytope)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
dim(::LazySet)
```
```@meta
CurrentModule = LazySets.VPolytopeModule
```
```@docs
dim(::VPolytope)
polyhedron(::VPolytope)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
rand(::Type{LazySet})
```
```@meta
CurrentModule = LazySets.VPolytopeModule
```
```@docs
rand(::Type{VPolytope})
remove_redundant_vertices(::VPolytope)
tohrep(::VPolytope)
tovrep(::VPolytope)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
linear_map(::AbstractMatrix, ::LazySet)
```
```@meta
CurrentModule = LazySets.VPolytopeModule
```
```@docs
linear_map(::AbstractMatrix, ::VPolytope)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
∈(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets.VPolytopeModule
```
```@docs
∈(::AbstractVector{N}, ::VPolytope{N}) where {N}
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
σ(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets.VPolytopeModule
```
```@docs
σ(::AbstractVector, ::VPolytope)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
cartesian_product(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets.VPolytopeModule
```
```@docs
cartesian_product(::VPolytope, ::VPolytope)
convex_hull(::VPolytope, ::VPolytope)
```
```@meta
CurrentModule = LazySets
```
```@docs; canonical=false
intersection(::Union{VPolygon,VPolytope}, ::Union{VPolygon,VPolytope})
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
minkowski_sum(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets.VPolytopeModule
```
```@docs
minkowski_sum(::VPolytope, ::VPolytope)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:
* [`extrema`](@ref extrema(::LazySet))
* [`extrema`](@ref extrema(::LazySet, ::Int))
* [`high`](@ref high(::LazySet))
* [`high`](@ref high(::LazySet, ::Int))
* [`isoperationtype`](@ref isoperationtype(::Type{LazySet}))
* [`low`](@ref low(::LazySet))
* [`low`](@ref low(::LazySet, ::Int))
* [`vertices_list`](@ref vertices_list(::LazySet))
* [`permute`](@ref permute(::LazySet, ::AbstractVector{Int}))
* [`project`](@ref project(::LazySet, ::AbstractVector{Int}))
* [`reflect`](@ref reflect(::LazySet))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`scale!`](@ref scale!(::Real, ::LazySet))
* [`ρ`](@ref ρ(::AbstractVector, ::LazySet))
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
* [`translate!`](@ref translate!(::LazySet, ::AbstractVector))

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
* [`norm`](@ref norm(::LazySet, ::Real=Inf))
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
* [`ispolyhedral`](@ref ispolyhedral(::AbstractPolyhedron))
* [`ispolyhedraltype`](@ref ispolyhedraltype(::Type{AbstractPolyhedron}))
* [`isdisjoint`](@ref isdisjoint(::AbstractPolyhedron, ::AbstractPolyhedron))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isboundedtype`](@ref isboundedtype(::Type{AbstractPolytope}))
* [`isempty`](@ref isempty(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope, ::Bool=false))
* [`volume`](@ref volume(::AbstractPolytope))
* [`⊆`](@ref ⊆(::AbstractPolytope, ::AbstractPolyhedron))
