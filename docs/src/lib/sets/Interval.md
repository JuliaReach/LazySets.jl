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
high(::Interval)
∈(::AbstractVector, ::Interval)
isflat(::Interval)
linear_map(::AbstractMatrix, ::Interval)
low(::Interval)
min(::Interval)
max(::Interval)
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
-(::Interval, ::Interval)
*(::Interval, ::Interval)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:
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
* [`generators`](@ref generators(::AbstractHyperrectangle))
* [`genmat`](@ref genmat(::AbstractHyperrectangle))
* [`high`](@ref high(::AbstractHyperrectangle, ::Int))
* [`low`](@ref low(::AbstractHyperrectangle, ::Int))
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))

Some additional functionality is available for `IntervalArithmetic.Interval`s:

```@docs
fast_interval_pow(::IA.Interval, ::Int)
```
