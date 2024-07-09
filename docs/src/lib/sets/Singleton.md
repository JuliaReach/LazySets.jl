```@meta
CurrentModule = LazySets.SingletonModule
```

# [Singleton](@id def_Singleton)

```@docs
Singleton
```

## Operations

```@docs
element(::Singleton)
rand(::Type{Singleton})
rectify(S::Singleton)
singleton_list(::Singleton)
linear_map(::AbstractMatrix, ::Singleton)
permute(::Singleton, ::AbstractVector{Int})
project(::Singleton, ::AbstractVector{Int})
translate!(::Singleton, ::AbstractVector)
```

```@meta
CurrentModule = LazySets.API
```

* [`isoperationtype`](@ref isoperationtype(::Type{LazySet}))
* [`scale!`](@ref scale!(::Real, ::LazySet))

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
* [`surface`](@ref surface(::LazySet))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`sample`](@ref sample(::LazySet))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
* [`≈`](@ref ≈(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))

Inherited from [`ConvexSet`](@ref):
* [`linear_combination`](@ref linear_combination(::ConvexSet, ::ConvexSet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`is_polyhedral`](@ref is_polyhedral(::AbstractPolyhedron))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isboundedtype`](@ref isboundedtype(::Type{AbstractPolytope}))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope))
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractCentrallySymmetricPolytope))

Inherited from [`AbstractZonotope`](@ref):
* [`order`](@ref order(::AbstractZonotope))
* [`togrep`](@ref togrep(::AbstractZonotope))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`area`](@ref area(::AbstractHyperrectangle))
* [`extrema`](@ref extrema(::AbstractHyperrectangle))
* [`extrema`](@ref extrema(::AbstractHyperrectangle, ::Int))
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real))
* [`volume`](@ref volume(::AbstractHyperrectangle))
* [`difference`](@ref difference(::AbstractHyperrectangle, ::AbstractHyperrectangle))

Inherited from [`AbstractSingleton`](@ref):
* [`center`](@ref center(::AbstractSingleton))
* [`center`](@ref center(::AbstractSingleton, ::Int))
* [`constraints_list`](@ref constraints_list(::AbstractSingleton))
* [`generators`](@ref generators(::AbstractSingleton{N}) where {N})
* [`genmat`](@ref genmat(::AbstractSingleton{N}) where {N})
* [`high`](@ref high(::AbstractSingleton))
* [`high`](@ref high(::AbstractSingleton, ::Int))
* [`isboundedtype`](@ref isboundedtype(::Type{AbstractSingleton}))
* [`low`](@ref low(::AbstractSingleton))
* [`low`](@ref low(::AbstractSingleton, ::Int))
* [`ngens`](@ref ngens(::AbstractSingleton))
* [`radius_hyperrectangle`](@ref radius_hyperrectangle(::AbstractSingleton{N}) where {N})
* [`radius_hyperrectangle`](@ref radius_hyperrectangle(::AbstractSingleton{N}, ::Int) where {N})
* [`reflect`](@ref reflect(::AbstractSingleton))
* [`vertices`](@ref vertices(::AbstractSingleton{N}) where {N})
* [`vertices_list`](@ref vertices_list(::AbstractSingleton))
* [`∈`](@ref ∈(::AbstractVector, ::AbstractSingleton))
* [`ρ`](@ref σ(::AbstractVector, ::AbstractSingleton))
* [`σ`](@ref σ(::AbstractVector, ::AbstractSingleton))
* [`cartesian_product`](@ref cartesian_product(::AbstractSingleton, ::AbstractSingleton))
* [`distance`](@ref distance(::AbstractSingleton, ::AbstractSingleton))
* [`intersection`](@ref intersection(::AbstractSingleton, ::AbstractSingleton))
* [`isdisjoint`](@ref isdisjoint(::AbstractSingleton, ::AbstractSingleton))
* [`isequivalent`](@ref isequivalent(::AbstractSingleton, ::AbstractSingleton))
* [`⊆`](@ref ⊆(::AbstractSingleton, ::AbstractSingleton))
* [`minkowski_difference`](@ref minkowski_difference(::AbstractSingleton, ::AbstractSingleton))
* [`minkowski_sum`](@ref minkowski_sum(::AbstractSingleton, ::AbstractSingleton))
