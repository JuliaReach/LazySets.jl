```@meta
CurrentModule = LazySets.IntervalModule
```

# [Interval](@id def_Interval)

```@docs
Interval
```

## Conversion

```@docs
convert(::Type{Interval}, ::LazySet)
convert(::Type{IA.Interval}, ::LazySet)
convert(::Type{Interval}, ::IA.Interval)
```

## Operations

```@docs
an_element(::Interval)
center(::Interval)
center(::Interval, ::Int)
chebyshev_center_radius(::Interval)
constraints_list(::Interval)
diameter(::Interval, ::Real=Inf)
dim(::Interval)
∈(::AbstractVector, ::Interval)
isflat(::Interval)
linear_map(::AbstractMatrix, ::Interval)
ngens(::Interval)
radius_hyperrectangle(::Interval)
radius_hyperrectangle(::Interval{N}, ::Int) where {N}
rand(::Type{Interval})
rectify(::Interval{N}) where {N}
reflect(::Interval)
scale(::Real, ::Interval)
split(::Interval, ::AbstractVector{Int})
ρ(::AbstractVector, ::Interval)
σ(::AbstractVector, ::Interval)
translate(::Interval, ::AbstractVector)
vertices_list(::Interval)
difference(::Interval{N}, ::Interval) where {N}
intersection(::Interval, ::Interval)
isdisjoint(::Interval, ::Interval, ::Bool=false)
⊆(::Interval, ::Interval, ::Bool=false)
minkowski_difference(::Interval, ::Interval)
minkowski_sum(::Interval, ::Interval)
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
* [`complement`](@ref complement(::LazySet))
* [`extrema`](@ref extrema(::LazySet))
* [`extrema`](@ref extrema(::LazySet, ::Int))
* [`high`](@ref high(::LazySet, ::Int))
* [`high`](@ref high(::LazySet))
* [`isoperationtype`](@ref isoperationtype(::Type{<:LazySet}))
* [`low`](@ref low(::LazySet, ::Int))
* [`low`](@ref low(::LazySet))
* [`norm`](@ref norm(::LazySet, ::Real))
* [`permute`](@ref permute(::LazySet, ::AbstractVector{Int}))
* [`project`](@ref project(::LazySet, ::AbstractVector{Int}))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`volume`](@ref volume(::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet, ::LazySet))
* [`≈`](@ref ≈(::LazySet, ::LazySet))

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
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`sample`](@ref sample(::LazySet, ::Int))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))

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
* [`distance`](@ref distance(::AbstractHyperrectangle, ::AbstractHyperrectangle))

Some additional functionality is available for `IntervalArithmetic.Interval`s:

```@docs
fast_interval_pow(::IA.Interval, ::Int)
```
