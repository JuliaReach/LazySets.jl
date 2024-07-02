```@meta
CurrentModule = LazySets.BallInfModule
```

# [Infinity-norm ball (BallInf)](@id def_BallInf)

```@docs
BallInf
```

## Operations

```@docs
center(::BallInf)
isflat(::BallInf)
ngens(::BallInf)
radius(::BallInf, ::Real=Inf)
radius_hyperrectangle(::BallInf)
radius_hyperrectangle(::BallInf, ::Int)
rand(::Type{BallInf})
reflect(::BallInf)
volume(::BallInf)
ρ(::AbstractVector, ::BallInf)
σ(::AbstractVector, ::BallInf)
translate!(::BallInf, ::AbstractVector)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:
* [`area`](@ref area(::LazySet))
* [`high`](@ref high(::LazySet))
* [`high`](@ref high(::LazySet, ::Int))
* [`isoperationtype`](@ref isoperationtype(::Type{LazySet}))
* [`low`](@ref low(::LazySet))
* [`low`](@ref low(::LazySet, ::Int))
* [`project`](@ref project(::LazySet, ::AbstractVector{Int}))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`minkowski_sum`](@ref minkowski_sum(::LazySet, ::LazySet))
```@meta
CurrentModule = LazySets
```
* [`ball_norm`](@ref ball_norm(::LazySet))
* [`radius_ball`](@ref radius_ball(::LazySet))

Inherited from [`LazySet`](@ref):
* [`complement`](@ref complement(::LazySet))
* [`concretize`](@ref concretize(::LazySet))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`constraints`](@ref constraints(::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet))
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`eltype`](@ref eltype(::LazySet))
* [`isoperation`](@ref isoperation(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`surface`](@ref surface(::LazySet))
* [`vertices`](@ref vertices(::LazySet))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
* [`sample`](@ref sample(::LazySet, ::Int))
* [`≈`](@ref ≈(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))

Inherited from [`LazySet`](@ref):
* [`linear_combination`](@ref linear_combination(::ConvexSet, ::ConvexSet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`is_polyhedral`](@ref is_polyhedral(::AbstractPolyhedron))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isboundedtype`](@ref isboundedtype(::Type{AbstractPolytope}))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope))
* [`center`](@ref center(::LazySet, ::Int))
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N})

Inherited from [`AbstractZonotope`](@ref):
* [`order`](@ref order(::AbstractZonotope))
* [`togrep`](@ref togrep(::AbstractZonotope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractZonotope))

Inherited from [`AbstractHyperrectangle`](@ref):
```@meta
CurrentModule = LazySets.API
```
* [`cartesian_product`](@ref cartesian_product(::LazySet, ::LazySet))
```@meta
CurrentModule = LazySets
```
* [`constraints_list`](@ref constraints_list(::AbstractHyperrectangle{N}) where {N})
* [`extrema`](@ref extrema(::AbstractHyperrectangle))
* [`extrema`](@ref extrema(::AbstractHyperrectangle, ::Int))
* [`generators`](@ref generators(::AbstractHyperrectangle))
* [`genmat`](@ref genmat(::AbstractHyperrectangle))
* [`isboundedtype`](@ref isboundedtype(::Type{AbstractHyperrectangle}))
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`rectify`](@ref rectify(::AbstractHyperrectangle))
* [`vertices_list`](@ref vertices_list(::AbstractHyperrectangle))
* [`∈`](@ref ∈(::AbstractVector, ::AbstractHyperrectangle))
* [`difference`](@ref difference(::AbstractHyperrectangle, ::AbstractHyperrectangle))
* [`distance`](@ref distance(::AbstractHyperrectangle, ::AbstractHyperrectangle))
* [`intersection`](@ref intersection(::AbstractHyperrectangle, ::AbstractHyperrectangle))
* [`isdisjoint`](@ref isdisjoint(::AbstractHyperrectangle, ::AbstractHyperrectangle))
* [`⊆`](@ref ⊆(::AbstractHyperrectangle, ::AbstractHyperrectangle))
* [`minkowski_difference`](@ref minkowski_difference(::AbstractHyperrectangle, ::AbstractHyperrectangle))
