```@meta
CurrentModule = LazySets.ZonotopeModule
```

# [Zonotope](@id def_Zonotope)

```@docs
Zonotope
```

## Conversion

```julia
convert(::Type{Zonotope}, ::AbstractZonotope)
```

## Operations

```@docs
generators(::Zonotope)
genmat(::Zonotope)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
rand(::Type{LazySet})
```
```@meta
CurrentModule = LazySets.ZonotopeModule
```
```@docs
rand(::Type{Zonotope})
```
```@meta
CurrentModule = LazySets
```
```@docs; canonical=false
remove_redundant_generators(::AbstractZonotope)
```
```@meta
CurrentModule = LazySets.ZonotopeModule
```
```@docs
remove_redundant_generators(::Zonotope)
remove_zero_generators(::Zonotope)
linear_map!(::Zonotope, ::AbstractMatrix, ::Zonotope)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
scale!(::Real, ::LazySet)
```
```@meta
CurrentModule = LazySets.ZonotopeModule
```
```@docs
scale!(::Real, ::Zonotope)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:

* [`center`](@ref center(::LazySet))
* [`high`](@ref high(::LazySet, ::Int))
* [`low`](@ref low(::LazySet, ::Int))
```@meta
CurrentModule = LazySets
```
* [`ngens`](@ref ngens(::AbstractZonotope))
```@meta
CurrentModule = LazySets.API
```
* [`isoperationtype`](@ref isoperationtype(::Type{LazySet}))
* [`permute`](@ref permute(::LazySet, ::AbstractVector{Int}))
```@meta
CurrentModule = LazySets
```
* [`reduce_order`](@ref reduce_order(::AbstractZonotope, ::Real, ::AbstractReductionMethod=GIR05()))
```@meta
CurrentModule = LazySets.API
```
* [`scale`](@ref scale(::Real, ::LazySet))
```@meta
CurrentModule = LazySets
```
* [`split`](@ref split(::AbstractZonotope, ::Int))
* [`split`](@ref split(::AbstractZonotope, ::AbstractVector{Int}, ::AbstractVector{Int}))
```@meta
CurrentModule = LazySets.API
```
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
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`eltype`](@ref eltype(::LazySet))
* [`isconvex`](@ref isconvex(::LazySet))
* [`isoperation`](@ref isoperation(::LazySet))
* [`ispolytopic`](@ref ispolytopic(::LazySet))
* [`ispolytopictype`](@ref ispolytopictype(::Type{LazySet}))
* [`norm`](@ref norm(::LazySet, ::Real))
* [`polyhedron`](@ref polyhedron(::LazySet))
* [`radius`](@ref radius(::LazySet, ::Real))
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
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
* [`convex_hull`](@ref convex_hull(::LazySet, ::LazySet))
* [`exact_sum`](@ref exact_sum(::LazySet, ::LazySet))
* [`isapprox`](@ref isapprox(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))

Inherited from [`ConvexSet`](@ref):
* [`isconvextype`](@ref isconvextype(::Type{ConvexSet}))
* [`linear_combination`](@ref linear_combination(::ConvexSet, ::ConvexSet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`ispolyhedral`](@ref ispolyhedral(::AbstractPolyhedron))
* [`ispolyhedraltype`](@ref ispolyhedraltype(::Type{AbstractPolyhedron}))
* [`intersection`](@ref intersection(::AbstractPolyhedron, ::LazySet))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`volume`](@ref volume(::AbstractPolytope))
* [`isboundedtype`](@ref isboundedtype(::Type{AbstractPolytope}))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope))
* [`center`](@ref center(::AbstractCentrallySymmetric, ::Int))
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetricPolytope))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetricPolytope, ::Int))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractCentrallySymmetricPolytope))

Inherited from [`AbstractZonotope`](@ref):
* [`constraints_list`](@ref constraints_list(::AbstractZonotope))
* [`constraints_list`](@ref constraints_list(::AbstractZonotope{<:AbstractFloat}; ::Bool=true))
* [`high`](@ref high(::AbstractZonotope))
* [`low`](@ref low(::AbstractZonotope))
* [`order`](@ref order(::AbstractZonotope))
* [`reflect`](@ref reflect(::AbstractZonotope))
* [`togrep`](@ref togrep(::AbstractZonotope))
* [`vertices_list`](@ref vertices_list(::AbstractZonotope))
* [`in`](@ref in(::AbstractVector, ::AbstractZonotope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractZonotope))
* [`project`](@ref project(::AbstractZonotope, ::AbstractVector{Int}))
* [`ρ`](@ref ρ(::AbstractVector, ::AbstractZonotope))
* [`σ`](@ref σ(::AbstractVector, ::AbstractZonotope))
* [`cartesian_product`](@ref cartesian_product(::AbstractZonotope, ::AbstractZonotope))
* [`isdisjoint`](@ref isdisjoint(::AbstractZonotope, ::AbstractZonotope))
* [`issubset`](@ref issubset(::AbstractZonotope, ::AbstractZonotope))
* [`minkowski_difference`](@ref minkowski_difference(::AbstractZonotope, ::AbstractZonotope))
* [`minkowski_sum`](@ref minkowski_sum(::AbstractZonotope, ::AbstractZonotope))
