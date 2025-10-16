```@meta
CurrentModule = LazySets.BallInfModule
```

# [Infinity-norm ball (BallInf)](@id def_BallInf)

```@docs
BallInf
```

## Operations

```@docs
isflat(::BallInf)
ngens(::BallInf)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
radius(::LazySet, ::Real=Inf)
```
```@meta
CurrentModule = LazySets.BallInfModule
```
```@docs
radius(::BallInf, ::Real=Inf)
radius_hyperrectangle(::BallInf)
radius_hyperrectangle(::BallInf, ::Int)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
rand(::Type{LazySet})
```
```@meta
CurrentModule = LazySets.BallInfModule
```
```@docs
rand(::Type{BallInf})
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
reflect(::LazySet)
```
```@meta
CurrentModule = LazySets.BallInfModule
```
```@docs
reflect(::BallInf)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
volume(::LazySet)
```
```@meta
CurrentModule = LazySets.BallInfModule
```
```@docs
volume(::BallInf)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
ρ(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets.BallInfModule
```
```@docs
ρ(::AbstractVector, ::BallInf)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
translate!(::LazySet, ::AbstractVector)
```
```@meta
CurrentModule = LazySets.BallInfModule
```
```@docs
translate!(::BallInf, ::AbstractVector)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:
* [`area`](@ref area(::LazySet))
* [`center`](@ref center(::LazySet))
* [`high`](@ref high(::LazySet))
* [`high`](@ref high(::LazySet, ::Int))
* [`isoperationtype`](@ref isoperationtype(::Type{LazySet}))
* [`low`](@ref low(::LazySet))
* [`low`](@ref low(::LazySet, ::Int))
* [`project`](@ref project(::LazySet, ::AbstractVector{Int}))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`σ`](@ref σ(::AbstractVector, ::LazySet))
* [`minkowski_sum`](@ref minkowski_sum(::LazySet, ::LazySet))
```@meta
CurrentModule = LazySets
```
* [`ball_norm`](@ref ball_norm(::LazySet))
* [`radius_ball`](@ref radius_ball(::LazySet))

Inherited from [`LazySet`](@ref):
* [`complement`](@ref complement(::LazySet))
* [`concretize`](@ref concretize(::LazySet))
* `copy(::Type{LazySet})`
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`constraints`](@ref constraints(::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet))
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`eltype`](@ref eltype(::LazySet))
* [`isoperation`](@ref isoperation(::LazySet))
* [`ispolytopic`](@ref ispolytopic(::LazySet))
* [`polyhedron`](@ref polyhedron(::LazySet))
* [`rationalize`](@ref rationalize(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`tosimplehrep`](@ref tosimplehrep(::LazySet))
* [`triangulate`](@ref triangulate(::LazySet))
* [`triangulate_faces`](@ref triangulate_faces(::LazySet))
* [`vertices`](@ref vertices(::LazySet))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`sample`](@ref sample(::LazySet, ::Int))
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
* [`convex_hull`](@ref convex_hull(::LazySet, ::LazySet))
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

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isboundedtype`](@ref isboundedtype(::Type{AbstractPolytope}))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope))
* [`center`](@ref center(::AbstractCentrallySymmetricPolytope, ::Int))
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
* [`constraints_list`](@ref constraints_list(::AbstractHyperrectangle))
* [`chebyshev_center_radius`](@ref chebyshev_center_radius(::AbstractHyperrectangle))
* [`distance`](@ref distance(::AbstractVector, ::AbstractHyperrectangle))
* [`extrema`](@ref extrema(::AbstractHyperrectangle))
* [`extrema`](@ref extrema(::AbstractHyperrectangle, ::Int))
* [`generators`](@ref generators(::AbstractHyperrectangle))
* [`genmat`](@ref genmat(::AbstractHyperrectangle))
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`rectify`](@ref rectify(::AbstractHyperrectangle))
* [`vertices_list`](@ref vertices_list(::AbstractHyperrectangle))
* [`∈`](@ref ∈(::AbstractVector, ::AbstractHyperrectangle))
* [`cartesian_product`](@ref cartesian_product(::AbstractHyperrectangle, ::AbstractHyperrectangle))
* [`difference`](@ref difference(::AbstractHyperrectangle, ::AbstractHyperrectangle))
* [`distance`](@ref distance(::AbstractHyperrectangle, ::AbstractHyperrectangle))
* [`intersection`](@ref intersection(::AbstractHyperrectangle, ::AbstractHyperrectangle))
* [`isdisjoint`](@ref isdisjoint(::AbstractHyperrectangle, ::AbstractHyperrectangle))
* [`⊆`](@ref ⊆(::AbstractHyperrectangle, ::AbstractHyperrectangle))
* [`minkowski_difference`](@ref minkowski_difference(::AbstractHyperrectangle, ::AbstractHyperrectangle))
