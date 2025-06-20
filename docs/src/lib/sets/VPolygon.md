```@meta
CurrentModule = LazySets.VPolygonModule
```

# [Polygon in vertex representation (VPolygon)](@id def_VPolygon)

```@docs
VPolygon
```

## Conversion

```julia
convert(::Type{VPolygon}, ::LazySet)
```

## Operations

```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
an_element(::LazySet)
```
```@meta
CurrentModule = LazySets.VPolygonModule
```
```@docs
an_element(::VPolygon)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
area(::LazySet)
```
```@meta
CurrentModule = LazySets.VPolygonModule
```
```@docs
area(::VPolygon)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
rand(::Type{LazySet})
```
```@meta
CurrentModule = LazySets.VPolygonModule
```
```@docs
rand(::Type{VPolygon})
remove_redundant_vertices(::VPolygon; ::String="monotone_chain")
remove_redundant_vertices!(::VPolygon; ::String="monotone_chain")
tohrep(::VPolygon, ::Type{HPOLYGON}=HPolygon) where {HPOLYGON<:AbstractHPolygon}
tovrep(::VPolygon)
∈(::AbstractVector, ::VPolygon)
linear_map(::AbstractMatrix, ::VPolygon)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
σ(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets.VPolygonModule
```
```@docs
σ(::AbstractVector, ::VPolygon)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
intersection(::LazySet, ::LazySet; ::Bool=true)
```
```@meta
CurrentModule = LazySets.VPolygonModule
```
```@docs
intersection(::VPolygon, ::VPolygon; ::Bool=true)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
minkowski_sum(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets.VPolygonModule
```
```@docs
minkowski_sum(::VPolygon, ::VPolygon)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:

* [`constraints_list`](@ref constraints_list(::LazySet))
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
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
* [`translate!`](@ref translate!(::LazySet, ::AbstractVector))
* [`convex_hull`](@ref convex_hull(::LazySet, ::LazySet))

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
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
* [`reflect`](@ref reflect(::LazySet))
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
* [`ρ`](@ref ρ(::AbstractVector, ::LazySet))
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
* [`ispolyhedral`](@ref ispolyhedral(::AbstractPolyhedron))
* [`ispolyhedraltype`](@ref ispolyhedraltype(::Type{AbstractPolyhedron}))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isboundedtype`](@ref isboundedtype(::Type{LazySet}))
* [`isempty`](@ref isempty(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope, ::Bool=false))
* [`isdisjoint`](@ref isdisjoint(::AbstractPolytope, ::AbstractPolytope))
* [`⊆`](@ref ⊆(::AbstractPolytope, ::AbstractPolyhedron))

Inherited from [`AbstractPolygon`](@ref):
* [`dim`](@ref dim(::AbstractPolygon))
* [`volume`](@ref volume(::AbstractPolygon))
