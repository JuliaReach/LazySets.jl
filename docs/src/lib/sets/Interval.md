```@meta
CurrentModule = LazySets
```

# [Interval](@id def_Interval)

```@docs
Interval
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
constraints_list(::Interval{N}) where {N}
rectify(::Interval{N}) where {N}
diameter(::Interval, ::Real=Inf)
split(::Interval, ::AbstractVector{Int})
chebyshev_center_radius(::Interval)
```

Inherited from [`ConvexSet`](@ref):
* [`singleton_list`](@ref singleton_list(::ConvexSet))

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
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real))
* [`ngens`](@ref ngens(::AbstractHyperrectangle))
* [`generators`](@ref generators(::AbstractHyperrectangle))
* [`genmat`](@ref genmat(::AbstractHyperrectangle))
* [`low`](@ref low(::AbstractHyperrectangle, ::Int))
* [`high`](@ref high(::AbstractHyperrectangle, ::Int))

Some additional functionality is available for `IntervalArithmetic.Interval`s:

```@docs
fast_interval_pow(::IA.Interval, ::Int)
```
