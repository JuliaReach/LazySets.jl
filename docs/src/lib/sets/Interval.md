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
dim(::Interval)
σ(::AbstractVector, ::Interval)
ρ(::AbstractVector, ::Interval)
∈(::AbstractVector, ::Interval)
∈(::Number, ::Interval)
an_element(::Interval)
vertices_list(::Interval)
translate(::Interval, ::AbstractVector)
center(::Interval)
center(::Interval, ::Int)
min(::Interval)
max(::Interval)
low(::Interval)
high(::Interval)
radius_hyperrectangle(::Interval)
radius_hyperrectangle(::Interval{N}, ::Int) where {N}
-(::Interval, ::Interval)
*(::Interval, ::Interval)
rand(::Type{Interval})
isflat(::Interval)
plot_recipe(::Interval{N}, ::Any=zero(N)) where {N}
linear_map(::AbstractMatrix, ::Interval)
scale(::Real, ::Interval)
constraints_list(::Interval)
rectify(::Interval{N}) where {N}
diameter(::Interval, ::Real=Inf)
split(::Interval, ::AbstractVector{Int})
ngens(::Interval)
chebyshev_center_radius(::Interval)
reflect(::Interval)
```

## Binary operations

```@docs
difference(::Interval{N}, ::Interval) where {N}
intersection(::Interval, ::Interval)
isdisjoint(::Interval, ::Interval, ::Bool=false)
⊆(::Interval, ::Interval, ::Bool=false)
minkowski_difference(::Interval, ::Interval)
minkowski_sum(::Interval, ::Interval)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:
* [`is_interior_point`](@ref is_interior_point(::AbstractVector{N}, ::LazySet; p=Inf, ε=_rtol(N)) where {N<:Real})
* [`radius`](@ref radius(::LazySet, ::Real))

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`singleton_list`](@ref singleton_list(::LazySet))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N})

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))

Inherited from [`AbstractZonotope`](@ref):
* [`order`](@ref order(::AbstractZonotope))
* [`togrep`](@ref togrep(::AbstractZonotope))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`generators`](@ref generators(::AbstractHyperrectangle))
* [`genmat`](@ref genmat(::AbstractHyperrectangle))
* [`low`](@ref low(::AbstractHyperrectangle, ::Int))
* [`high`](@ref high(::AbstractHyperrectangle, ::Int))

Some additional functionality is available for `IntervalArithmetic.Interval`s:

```@docs
fast_interval_pow(::IA.Interval, ::Int)
```
