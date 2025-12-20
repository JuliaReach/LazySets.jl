```@meta
CurrentModule = LazySets.IntervalModule
```

# [Interval](@id def_Interval)

```@docs
Interval
```

## Conversion

```julia
convert(::Type{Interval}, ::LazySet)
convert(::Type{IA.Interval}, ::LazySet)
convert(::Type{Interval}, ::IA.Interval)
```

## Operations

```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
an_element(::LazySet)
```
```@meta
CurrentModule = LazySets.IntervalModule
```
```@docs
an_element(::Interval)
chebyshev_center_radius(::Interval)
isflat(::Interval)
ngens(::Interval)
radius_hyperrectangle(::Interval)
radius_hyperrectangle(::Interval, ::Int)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
rand(::Type{<:LazySet})
```
```@meta
CurrentModule = LazySets.IntervalModule
```
```@docs
rand(::Type{Interval})
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
rectify(::LazySet)
```
```@meta
CurrentModule = LazySets.IntervalModule
```
```@docs
rectify(::Interval)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
vertices_list(::LazySet)
```
```@meta
CurrentModule = LazySets.IntervalModule
```
```@docs
vertices_list(::Interval)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
linear_map(::AbstractMatrix, ::LazySet)
```
```@meta
CurrentModule = LazySets.IntervalModule
```
```@docs
linear_map(::AbstractMatrix, ::Interval)
split(::Interval, ::AbstractVector{Int})
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
translate(::LazySet, ::AbstractVector)
```
```@meta
CurrentModule = LazySets.IntervalModule
```
```@docs
translate(::Interval, ::AbstractVector)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
difference(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets.IntervalModule
```
```@docs
difference(::Interval, ::Interval)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
isdisjoint(::LazySet, ::LazySet, ::Bool=false)
```
```@meta
CurrentModule = LazySets.IntervalModule
```
```@docs
isdisjoint(::Interval, ::Interval, ::Bool=false)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
minkowski_difference(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets.IntervalModule
```
```@docs
minkowski_difference(::Interval, ::Interval)
plot_recipe(::Interval{N}, ::Any=zero(N)) where {N}
min(::Interval)
max(::Interval)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:
* [`center`](@ref center(::LazySet))
* [`center`](@ref center(::LazySet, ::Int))
* [`complement`](@ref complement(::LazySet))
* [`constraints_list`](@ref constraints_list(::LazySet))
* `copy(::Type{LazySet})`
* [`diameter`](@ref diameter(::LazySet, ::Real=Inf))
* [`dim`](@ref dim(::LazySet))
* [`extrema`](@ref extrema(::LazySet))
* [`extrema`](@ref extrema(::LazySet, ::Int))
* [`high`](@ref high(::LazySet, ::Int))
* [`high`](@ref high(::LazySet))
* [`isoperationtype`](@ref isoperationtype(::Type{<:LazySet}))
* [`low`](@ref low(::LazySet, ::Int))
* [`low`](@ref low(::LazySet))
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`reflect`](@ref reflect(::LazySet))
* [`volume`](@ref volume(::LazySet))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`in`](@ref in(::AbstractVector, ::LazySet))
* [`permute`](@ref permute(::LazySet, ::AbstractVector{Int}))
* [`project`](@ref project(::LazySet, ::AbstractVector{Int}))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`ρ`](@ref ρ(::AbstractVector, ::LazySet))
* [`σ`](@ref σ(::AbstractVector, ::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet, ::LazySet))
* [`distance`](@ref distance(::LazySet, ::LazySet))
* [`intersection`](@ref intersection(::LazySet, ::LazySet))
* [`isapprox`](@ref isapprox(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))
* [`issubset`](@ref issubset(::LazySet, ::LazySet, ::Bool=false))
* [`minkowski_sum`](@ref minkowski_sum(::LazySet, ::LazySet))

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`concretize`](@ref concretize(::LazySet))
* [`constraints`](@ref constraints(::LazySet))
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`eltype`](@ref eltype(::LazySet))
* [`isconvex`](@ref isconvex(::LazySet))
* [`isoperation`](@ref isoperation(::LazySet))
* [`ispolytopic`](@ref ispolytopic(::LazySet))
* [`ispolytopictype`](@ref ispolytopictype(::Type{LazySet}))
* [`polyhedron`](@ref polyhedron(::LazySet))
* [`rationalize`](@ref rationalize(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`tosimplehrep`](@ref tosimplehrep(::LazySet))
* [`vertices`](@ref vertices(::LazySet))
* [`sample`](@ref sample(::LazySet, ::Int))
* [`==`](@ref ==(::LazySet, ::LazySet))

Inherited from [`LazySet`](@ref) but does not apply:
* [`triangulate`](@ref triangulate(::LazySet))
* [`triangulate_faces`](@ref triangulate_faces(::LazySet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`ispolyhedral`](@ref ispolyhedral(::AbstractPolyhedron))
* [`ispolyhedraltype`](@ref ispolyhedraltype(::Type{AbstractPolyhedron}))

Inherited from [`ConvexSet`](@ref):
* [`isconvextype`](@ref isconvextype(::Type{ConvexSet}))
* [`linear_combination`](@ref linear_combination(::ConvexSet, ::ConvexSet))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isboundedtype`](@ref isboundedtype(::Type{<:AbstractPolytope}))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractCentrallySymmetricPolytope, ::Bool=false))

Inherited from [`AbstractZonotope`](@ref):
* [`order`](@ref order(::AbstractZonotope))
* [`togrep`](@ref togrep(::AbstractZonotope))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`generators`](@ref generators(::AbstractHyperrectangle))
* [`genmat`](@ref genmat(::AbstractHyperrectangle))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::AbstractHyperrectangle))
* [`cartesian_product`](@ref cartesian_product(::AbstractHyperrectangle, ::AbstractHyperrectangle))

Some additional functionality is available for `IntervalArithmetic.Interval`s:

```@docs
fast_interval_pow(::IA.Interval, ::Int)
```
