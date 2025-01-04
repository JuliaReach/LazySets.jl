```@contents
Pages = ["AbstractSingleton.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# [Singletons (AbstractSingleton)](@id def_AbstractSingleton)

A singleton is a special hyperrectangle consisting of only one point.

```@docs
AbstractSingleton
```

This interface requires to implement the following function:

```@docs
element(::AbstractSingleton)
```

This interface defines the following functions:

```@docs
element(::AbstractSingleton, ::Int)
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
reflect(::AbstractSingleton)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
∈(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
∈(::AbstractVector, ::AbstractSingleton)
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
σ(::AbstractVector, ::AbstractSingleton)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:

* [`center`](@ref center(::LazySet))
* [`center`](@ref center(::LazySet, ::Int))
```@meta
CurrentModule = LazySets
```
* [`chebyshev_center_radius`](@ref chebyshev_center_radius(::LazySet))
```@meta
CurrentModule = LazySets.API
```
* [`constraints_list`](@ref constraints_list(::LazySet))
```@meta
CurrentModule = LazySets
```
* [`generators`](@ref generators(::AbstractZonotope))
* [`genmat`](@ref genmat(::AbstractZonotope))
```@meta
CurrentModule = LazySets.API
```
* [`high`](@ref high(::LazySet))
* [`high`](@ref high(::LazySet, ::Int))
* [`low`](@ref low(::LazySet))
* [`low`](@ref low(::LazySet, ::Int))
```@meta
CurrentModule = LazySets
```
* [`ngens`](@ref ngens(::AbstractZonotope))
* [`radius_hyperrectangle`](@ref radius_hyperrectangle(::AbstractVector, ::AbstractHyperrectangle))
* [`radius_hyperrectangle`](@ref radius_hyperrectangle(::AbstractVector, ::AbstractHyperrectangle, ::Int))
```@meta
CurrentModule = LazySets.API
```
* [`vertices`](@ref vertices(::LazySet))
* [`vertices_list`](@ref vertices_list(::LazySet))
* [`ρ`](@ref ρ(::AbstractVector, ::LazySet))
* [`cartesian_product`](@ref cartesian_product(::LazySet, ::LazySet))
* [`distance`](@ref distance(::LazySet, ::LazySet; ::Real=2.0))
* [`intersection`](@ref intersection(::LazySet, ::LazySet))
* [`isdisjoint`](@ref isdisjoint(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊆`](@ref ⊆(::LazySet, ::LazySet))
* [`minkowski_sum`](@ref minkowski_sum(::LazySet, ::LazySet))

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
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
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`sample`](@ref sample(::LazySet, ::Int=1))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
* [`convex_hull`](@ref convex_hull(::LazySet, ::LazySet))
* [`exact_sum`](@ref exact_sum(::LazySet, ::LazySet))
* [`≈`](@ref ≈(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))

Inherited from [`ConvexSet`](@ref):
* [`isconvextype`](@ref isconvextype(::Type{LazySet}))
* [`linear_combination`](@ref linear_combination(::ConvexSet, ::ConvexSet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`ispolyhedral`](@ref ispolyhedral(::AbstractPolyhedron))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isboundedtype`](@ref isboundedtype(::Type{AbstractPolytope}))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope))
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractCentrallySymmetricPolytope, ::Bool=false))

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
* [`extrema`](@ref extrema(::AbstractHyperrectangle))
* [`extrema`](@ref extrema(::AbstractHyperrectangle, ::Int))
* [`isflat`](@ref isflat(::AbstractHyperrectangle))
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real=Inf))
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real=Inf))
* [`rectify`](@ref rectify(::AbstractHyperrectangle))
* [`volume`](@ref volume(::AbstractHyperrectangle))
* [`distance`](@ref distance(::AbstractVector, ::AbstractHyperrectangle; ::Real=2.0))
* [`difference`](@ref difference(::AbstractHyperrectangle, ::AbstractHyperrectangle))
* [`minkowski_difference`](@ref minkowski_difference(::AbstractHyperrectangle, ::AbstractHyperrectangle))

Plotting singletons is available too:

```@docs
plot_recipe(::AbstractSingleton{N}, ::Any=zero(N)) where {N}
RecipesBase.apply_recipe(::AbstractDict{Symbol,Any}, ::AbstractSingleton{N}, ::Real=zero(N)) where {N}
```

## Implementations

* [Singleton](@ref def_Singleton)
* [Origin (ZeroSet)](@ref def_ZeroSet)
