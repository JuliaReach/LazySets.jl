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

```@docs; canonical=false
an_element(::LazySet)
```
```@docs
an_element(::Interval)
chebyshev_center_radius(::Interval)
isflat(::Interval)
```
```@docs; canonical=false
linear_map(::AbstractMatrix, ::LazySet)
```
```@docs
linear_map(::AbstractMatrix, ::Interval)
ngens(::Interval)
radius_hyperrectangle(::Interval)
radius_hyperrectangle(::Interval{N}, ::Int) where {N}
```
```@docs; canonical=false
rand(::Type{LazySet})
```
```@docs
rand(::Type{Interval})
```
```@docs; canonical=false
rectify(::LazySet, ::AbstractVector{Int})
```
```@docs
rectify(::Interval{N}) where {N}
split(::Interval, ::AbstractVector{Int})
```
```@docs; canonical=false
translate(::LazySet, ::AbstractVector)
```
```@docs
translate(::Interval, ::AbstractVector)
```
```@docs; canonical=false
vertices_list(::LazySet)
```
```@docs
vertices_list(::Interval)
```
```@docs; canonical=false
difference(::LazySet, ::LazySet)
```
```@docs
difference(::Interval{N}, ::Interval) where {N}
```
```@docs; canonical=false
isdisjoint(::LazySet, ::LazySet, ::Bool=false)
```
```@docs
isdisjoint(::Interval, ::Interval, ::Bool=false)
```
```@docs; canonical=false
minkowski_difference(::LazySet, ::LazySet)
```
```@docs
minkowski_difference(::Interval, ::Interval)
plot_recipe(::Interval{N}, ::Any=zero(N)) where {N}
min(::Interval)
max(::Interval)
-(::Interval, ::Interval)
*(::Interval, ::Interval)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`center`](@ref center(::LazySet))
* [`center`](@ref center(::LazySet, ::Int))
* [`complement`](@ref complement(::LazySet))
* [`constraints_list`](@ref constraints_list(::LazySet))
* [`diameter`](@ref diameter(::LazySet, ::Real=Inf))
* [`dim`](@ref dim(::LazySet))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`extrema`](@ref extrema(::LazySet))
* [`extrema`](@ref extrema(::LazySet, ::Int))
* [`high`](@ref high(::LazySet, ::Int))
* [`high`](@ref high(::LazySet))
* [`∈`](@ref ∈(::AbstractVector, ::LazySet))
* [`isoperationtype`](@ref isoperationtype(::Type{<:LazySet}))
* [`low`](@ref low(::LazySet, ::Int))
* [`low`](@ref low(::LazySet))
* [`norm`](@ref norm(::LazySet, ::Real))
* [`permute`](@ref permute(::LazySet, ::AbstractVector{Int}))
* [`project`](@ref project(::LazySet, ::AbstractVector{Int}))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`reflect`](@ref reflect(::LazySet))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`ρ`](@ref ρ(::AbstractVector, ::LazySet))
* [`σ`](@ref σ(::AbstractVector, ::LazySet))
* [`volume`](@ref volume(::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet, ::LazySet))
* [`distance`](@ref distance(::LazySet, ::LazySet))
* [`intersection`](@ref intersection(::LazySet, ::LazySet))
* [`≈`](@ref ≈(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))
* [`⊆`](@ref ⊆(::LazySet, ::LazySet, ::Bool=false))
* [`minkowski_sum`](@ref minkowski_sum(::LazySet, ::LazySet))

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`concretize`](@ref concretize(::LazySet))
* [`constraints`](@ref constraints(::LazySet))
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`eltype`](@ref eltype(::LazySet))
* [`isoperation`](@ref isoperation(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`vertices`](@ref vertices(::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`sample`](@ref sample(::LazySet, ::Int))
* [`==`](@ref ==(::LazySet, ::LazySet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`is_polyhedral`](@ref is_polyhedral(::AbstractPolyhedron))

Inherited from [`ConvexSet`](@ref):

```@meta
CurrentModule = LazySets.API
```

* [`linear_combination`](@ref linear_combination(::LazySet, ::LazySet))

```@meta
CurrentModule = LazySets
```

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isboundedtype`](@ref isboundedtype(::Type{<:AbstractPolytope}))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractCentrallySymmetricPolytope{N}, ::Bool=false) where {N})

Inherited from [`AbstractZonotope`](@ref):
* [`order`](@ref order(::AbstractZonotope))
* [`togrep`](@ref togrep(::AbstractZonotope))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`generators`](@ref generators(::AbstractHyperrectangle))
* [`genmat`](@ref genmat(::AbstractHyperrectangle))
* [`isconvextype`](@ref isconvextype(::Type{<:AbstractHyperrectangle}))
* [`cartesian_product`](@ref cartesian_product(::AbstractHyperrectangle, ::AbstractHyperrectangle))

Some additional functionality is available for `IntervalArithmetic.Interval`s:

```@docs
fast_interval_pow(::IA.Interval, ::Int)
```
