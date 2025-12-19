```@meta
CurrentModule = LazySets.StarModule
```

# [Star](@id def_Star)

```@docs
Star
```

## Conversion

```julia
convert(::Type{Star}, ::AbstractPolyhedron)
```

## Operations

```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
an_element(::LazySet)
```
```@meta
CurrentModule = LazySets.StarModule
```
```@docs
an_element(::Star)
basis(::Star)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
constraints_list(::LazySet)
```
```@meta
CurrentModule = LazySets.StarModule
```
```@docs
constraints_list(::Star)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
isbounded(::LazySet)
```
```@meta
CurrentModule = LazySets.StarModule
```
```@docs
isbounded(::Star)
predicate(::Star)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
rand(::Type{LazySet})
```
```@meta
CurrentModule = LazySets.StarModule
```
```@docs
rand(::Type{Star})
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
in(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets.StarModule
```
```@docs
in(::AbstractVector, ::Star)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:

* [`center`](@ref center(::LazySet))
* [`center`](@ref center(::LazySet, ::Int))
* [`dim`](@ref dim(::LazySet))
* [`isempty`](@ref isempty(::LazySet))
* [`isoperationtype`](@ref isoperationtype(::Type{LazySet}))
* [`vertices_list`](@ref vertices_list(::LazySet))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::LazySet))
* [`ρ`](@ref ρ(::AbstractVector, ::LazySet))
* [`σ`](@ref σ(::AbstractVector, ::LazySet))

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
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`sample`](@ref sample(::LazySet, ::Int=1))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
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
* [`isuniversal`](@ref isuniversal(::AbstractPolyhedron, ::Bool=false))
* [`low`](@ref low(::AbstractPolyhedron))
* [`low`](@ref low(::AbstractPolyhedron, ::Int))
* [`project`](@ref project(::AbstractPolyhedron, ::AbstractVector{Int}))
* [`intersection`](@ref intersection(::AbstractPolyhedron, ::AbstractPolyhedron))
* [`isdisjoint`](@ref isdisjoint(::AbstractPolyhedron, ::AbstractPolyhedron))
* [`issubset`](@ref issubset(::AbstractPolyhedron, ::AbstractPolyhedron))
* [`minkowski_sum`](@ref minkowski_sum(::AbstractPolyhedron, ::AbstractPolyhedron))
