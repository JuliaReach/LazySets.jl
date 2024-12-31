```@meta
CurrentModule = LazySets.ZeroSetModule
```

# [Origin (ZeroSet)](@id def_ZeroSet)

```@docs
ZeroSet
```

## Operations

```@docs
element(::ZeroSet)
element(::ZeroSet, ::Int)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
rand(::Type{LazySet})
```
```@meta
CurrentModule = LazySets.ZeroSetModule
```
```@docs
rand(::Type{ZeroSet})
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:

* `copy(::Type{ZeroSet})`
* [`dim`](@ref dim(::LazySet))
* [`reflect`](@ref reflect(::LazySet))
* [`rectify`](@ref rectify(::LazySet))
* [`∈`](@ref ∈(::AbstractVector, ::LazySet))
* [`isoperationtype`](@ref isoperationtype(::Type{<:LazySet}))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::LazySet))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`scale!`](@ref scale!(::Real, ::LazySet))
* [`ρ`](@ref ρ(::AbstractVector, ::LazySet))
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
* [`minkowski_sum`](@ref minkowski_sum(::LazySet, ::LazySet))

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`complement`](@ref complement(::LazySet))
* [`concretize`](@ref concretize(::LazySet))
* [`constraints`](@ref constraints(::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`eltype`](@ref eltype(::LazySet))
* [`isoperation`](@ref isoperation(::LazySet))
* [`rationalize`](@ref rationalize(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`surface`](@ref surface(::LazySet))
* [`triangulate`](@ref triangulate(::LazySet))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`sample`](@ref sample(::LazySet, ::Int))
* [`≈`](@ref ≈(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))

Inherited from [`ConvexSet`](@ref):
* [`isconvextype`](@ref isconvextype(::Type{ConvexSet}))
* [`linear_combination`](@ref linear_combination(::ConvexSet, ::ConvexSet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`ispolyhedral`](@ref ispolyhedral(::AbstractPolyhedron))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isboundedtype`](@ref isboundedtype(::Type{AbstractPolytope}))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractCentrallySymmetricPolytope, ::Bool=false))

Inherited from [`AbstractZonotope`](@ref):
* [`order`](@ref order(::AbstractZonotope))
* [`remove_redundant_generators`](@ref remove_redundant_generators(::AbstractZonotope))
* [`togrep`](@ref togrep(::AbstractZonotope))
* [`reduce_order`](@ref reduce_order(::AbstractZonotope, ::Real, ::AbstractReductionMethod=GIR05()))
* [`split`](@ref split(::AbstractZonotope, ::Int))
* [`split`](@ref split(::AbstractZonotope, ::AbstractVector{Int}, ::AbstractVector{Int}))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`area`](@ref area(::AbstractHyperrectangle))
* [`extrema`](@ref extrema(::AbstractHyperrectangle))
* [`extrema`](@ref extrema(::AbstractHyperrectangle, ::Int))
* [`isflat`](@ref isflat(::AbstractHyperrectangle))
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real))
* [`volume`](@ref volume(::AbstractHyperrectangle))
* [`project`](@ref project(::AbstractHyperrectangle, ::AbstractVector{Int}))
* [`difference`](@ref difference(::AbstractHyperrectangle, ::AbstractHyperrectangle))

Inherited from [`AbstractSingleton`](@ref):
* [`center`](@ref center(::AbstractSingleton))
* [`center`](@ref center(::AbstractSingleton, ::Int))
* [`chebyshev_center_radius`](@ref chebyshev_center_radius(::AbstractSingleton))
* [`constraints_list`](@ref constraints_list(::AbstractSingleton))
* [`generators`](@ref generators(::AbstractSingleton))
* [`genmat`](@ref genmat(::AbstractSingleton))
* [`high`](@ref high(::AbstractSingleton))
* [`high`](@ref high(::AbstractSingleton, ::Int))
* [`low`](@ref low(::AbstractSingleton))
* [`low`](@ref low(::AbstractSingleton, ::Int))
* [`ngens`](@ref ngens(::AbstractSingleton))
* [`radius_hyperrectangle`](@ref radius_hyperrectangle(::AbstractSingleton))
* [`radius_hyperrectangle`](@ref radius_hyperrectangle(::AbstractSingleton, ::Int))
* [`vertices`](@ref vertices(::AbstractSingleton))
* [`vertices`](@ref vertices(::AbstractSingleton))
* [`vertices_list`](@ref vertices_list(::AbstractSingleton))
* [`σ`](@ref σ(::AbstractVector, ::AbstractSingleton))
* [`cartesian_product`](@ref cartesian_product(::AbstractSingleton, ::AbstractSingleton))
* [`distance`](@ref distance(::AbstractSingleton, ::AbstractSingleton))
* [`intersection`](@ref intersection(::AbstractSingleton, ::AbstractSingleton))
* [`isdisjoint`](@ref isdisjoint(::AbstractSingleton, ::AbstractSingleton))
* [`isequivalent`](@ref isequivalent(::AbstractSingleton, ::AbstractSingleton))
* [`⊆`](@ref ⊆(::AbstractSingleton, ::AbstractSingleton))
