```@meta
CurrentModule = LazySets.HyperrectangleModule
```

# [Hyperrectangle](@id def_Hyperrectangle)

```@docs
Hyperrectangle
```

## Conversion

```julia
convert(::Type{Hyperrectangle}, ::AbstractHyperrectangle)
convert(::Type{Hyperrectangle}, ::IA.IntervalBox)
convert(::Type{IA.IntervalBox}, ::AbstractHyperrectangle)
```

## Operations

```@docs
radius_hyperrectangle(::Hyperrectangle)
radius_hyperrectangle(::Hyperrectangle, ::Int)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
rand(::Type{LazySet})
```
```@meta
CurrentModule = LazySets.HyperrectangleModule
```
```@docs
rand(::Type{Hyperrectangle})
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
translate(::LazySet, ::AbstractVector)
```
```@meta
CurrentModule = LazySets.HyperrectangleModule
```
```@docs
translate(::Hyperrectangle, ::AbstractVector)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:

* [`center`](@ref center(::LazySet))
```@meta
CurrentModule = LazySets
```
* [`genmat`](@ref genmat(::AbstractHyperrectangle))
```@meta
CurrentModule = LazySets.API
```
* [`isoperationtype`](@ref isoperationtype(::Type{LazySet}))
* [`permute`](@ref permute(::LazySet, ::AbstractVector{Int}))
* [`scale!`](@ref scale!(::Real, ::LazySet))
* [`ρ`](@ref ρ(::AbstractVector, ::LazySet))
* [`σ`](@ref σ(::AbstractVector, ::LazySet))
* [`translate!`](@ref translate!(::LazySet, ::AbstractVector))

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
* [`delaunay`](@ref delaunay(::LazySet))
* [`diameter`](@ref diameter(::LazySet, ::Real=Inf))
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`eltype`](@ref eltype(::LazySet))
* [`isoperation`](@ref isoperation(::LazySet))
* [`polyhedron`](@ref polyhedron(::LazySet))
* [`rationalize`](@ref rationalize(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`surface`](@ref surface(::LazySet))
* [`tosimplehrep`](@ref tosimplehrep(::LazySet))
* [`triangulate`](@ref triangulate(::LazySet))
* [`vertices`](@ref vertices(::LazySet))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`sample`](@ref sample(::LazySet, ::Int=1))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet, ::LazySet))
* [`exact_sum`](@ref exact_sum(::LazySet, ::LazySet))
* [`≈`](@ref ≈(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))

Inherited from [`ConvexSet`](@ref):
* [`isconvextype`](@ref isconvextype(::Type{ConvexSet}))
* [`linear_combination`](@ref linear_combination(::ConvexSet, ::ConvexSet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`ispolyhedral`](@ref ispolyhedral(::AbstractPolyhedron))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isboundedtype`](@ref isboundedtype(::Type{LazySet}))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope))
* [`center`](@ref center(::AbstractCentrallySymmetricPolytope, ::Int))
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractCentrallySymmetricPolytope))

Inherited from [`AbstractZonotope`](@ref):
* [`order`](@ref order(::AbstractZonotope))
* [`remove_redundant_generators`](@ref remove_redundant_generators(::AbstractZonotope))
* [`togrep`](@ref togrep(::AbstractZonotope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractZonotope))
* [`reduce_order`](@ref reduce_order(::AbstractZonotope, ::Real, ::AbstractReductionMethod=GIR05()))
* [`split`](@ref split(::AbstractZonotope, ::Int))
* [`split`](@ref split(::AbstractZonotope, ::AbstractVector{Int}, ::AbstractVector{Int}))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`area`](@ref area(::AbstractHyperrectangle))
* [`constraints_list`](@ref constraints_list(::AbstractHyperrectangle))
* [`extrema`](@ref extrema(::AbstractHyperrectangle))
* [`extrema`](@ref extrema(::AbstractHyperrectangle, ::Int))
* [`generators`](@ref generators(::AbstractHyperrectangle))
* [`genmat`](@ref genmat(::AbstractHyperrectangle))
* [`high`](@ref high(::AbstractHyperrectangle))
* [`high`](@ref high(::AbstractHyperrectangle, ::Int))
* [`isflat`](@ref isflat(::Hyperrectangle))
* [`low`](@ref low(::AbstractHyperrectangle))
* [`low`](@ref low(::AbstractHyperrectangle, ::Int))
* [`ngens`](@ref ngens(::AbstractHyperrectangle))
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real))
* [`rectify`](@ref rectify(::AbstractHyperrectangle))
* [`reflect`](@ref reflect(::AbstractHyperrectangle))
* [`vertices_list`](@ref vertices_list(::AbstractHyperrectangle))
* [`volume`](@ref volume(::AbstractHyperrectangle))
* [`distance`](@ref distance(::AbstractVector, ::AbstractHyperrectangle))
* [`∈`](@ref ∈(::AbstractVector, ::AbstractHyperrectangle))
* [`project`](@ref project(::AbstractHyperrectangle, ::AbstractVector{Int}))
* [`split`](@ref split(::AbstractHyperrectangle, ::AbstractVector{Int}))
* [`ρ`](@ref ρ(::AbstractVector, ::AbstractHyperrectangle))
* [`σ`](@ref σ(::AbstractVector, ::AbstractHyperrectangle))
* [`cartesian_product`](@ref cartesian_product(::AbstractHyperrectangle, ::AbstractHyperrectangle))
* [`difference`](@ref difference(::AbstractHyperrectangle, ::AbstractHyperrectangle))
* [`distance`](@ref distance(::AbstractHyperrectangle, ::AbstractHyperrectangle; ::Real=2.0))
* [`intersection`](@ref intersection(::AbstractHyperrectangle, ::AbstractHyperrectangle))
* [`isdisjoint`](@ref isdisjoint(::AbstractHyperrectangle, ::AbstractHyperrectangle))
* [`⊆`](@ref ⊆(::AbstractHyperrectangle, ::AbstractHyperrectangle))
* [`minkowski_difference`](@ref minkowski_difference(::AbstractHyperrectangle, ::AbstractHyperrectangle))
* [`minkowski_sum`](@ref minkowski_sum(::AbstractHyperrectangle, ::AbstractHyperrectangle))
