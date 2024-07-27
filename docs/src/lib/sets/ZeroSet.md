```@meta
CurrentModule = LazySets.ZeroSetModule
```

# [Origin (ZeroSet)](@id def_ZeroSet)

```@docs
ZeroSet
```

## Operations

```@docs
dim(::ZeroSet)
element(::ZeroSet{N}) where {N}
element(::ZeroSet{N}, ::Int) where {N}
rand(::Type{ZeroSet})
rectify(::ZeroSet)
reflect(::ZeroSet)
∈(::AbstractVector, ::ZeroSet)
linear_map(::AbstractMatrix, ::ZeroSet)
ρ(::AbstractVector, ::ZeroSet)
translate(::ZeroSet, ::AbstractVector)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:
* [`isoperationtype`](@ref isoperationtype(::Type{<:LazySet}))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`scale!`](@ref scale!(::Real, ::LazySet))
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
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`surface`](@ref surface(::LazySet))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`sample`](@ref sample(::LazySet, ::Int))
* [`≈`](@ref ≈(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))

Inherited from [`ConvexSet`](@ref):
* [`linear_combination`](@ref linear_combination(::ConvexSet, ::ConvexSet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`ispolyhedral`](@ref ispolyhedral(::AbstractPolyhedron))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isboundedtype`](@ref isboundedtype(::Type{AbstractPolytope}))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope))
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
* [`isconvextype`](@ref isconvextype(::Type{AbstractSingleton}))
* [`low`](@ref low(::AbstractSingleton))
* [`low`](@ref low(::AbstractSingleton, ::Int))
* [`ngens`](@ref ngens(::AbstractSingleton))
* [`radius_hyperrectangle`](@ref radius_hyperrectangle(::AbstractSingleton{N}) where {N})
* [`radius_hyperrectangle`](@ref radius_hyperrectangle(::AbstractSingleton{N}, ::Int) where {N})
* [`vertices`](@ref vertices(::AbstractSingleton{N}) where {N})
* [`vertices`](@ref vertices(::AbstractSingleton))
* [`vertices_list`](@ref vertices_list(::AbstractSingleton))
* [`σ`](@ref σ(::AbstractVector, ::AbstractSingleton))
* [`cartesian_product`](@ref cartesian_product(::AbstractSingleton, ::AbstractSingleton))
* [`distance`](@ref distance(::AbstractSingleton, ::AbstractSingleton))
* [`intersection`](@ref intersection(::AbstractSingleton, ::AbstractSingleton))
* [`isdisjoint`](@ref isdisjoint(::AbstractSingleton, ::AbstractSingleton))
* [`isequivalent`](@ref isequivalent(::AbstractSingleton, ::AbstractSingleton))
* [`⊆`](@ref ⊆(::AbstractSingleton, ::AbstractSingleton))
