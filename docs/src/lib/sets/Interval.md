```@meta
CurrentModule = LazySets
```

# [Interval](@id def_Interval)

```@docs
Interval
dim(::Interval)
σ(::AbstractVector{N}, ::Interval{N}) where {N<:Real}
ρ(::AbstractVector{N}, ::Interval{N}) where {N<:Real}
∈(::AbstractVector{N}, ::Interval{N}) where {N<:Real}
∈(::N, ::Interval{N}) where {N<:Real}
an_element(::Interval{N}) where {N<:Real}
vertices_list(::Interval{N}) where {N<:Real}
translate(::Interval{N}, ::AbstractVector{N}) where {N<:Real}
center(::Interval{N}) where {N<:Real}
min(::Interval{N}) where {N<:Real}
max(::Interval{N}) where {N<:Real}
low(::Interval{N}) where {N<:Real}
high(::Interval{N}) where {N<:Real}
radius_hyperrectangle(::Interval{N}) where {N<:Real}
radius_hyperrectangle(::Interval{N}, ::Int) where {N<:Real}
-(::Interval{N}, ::Interval{N}) where {N<:Real}
*(::Interval{N}, ::Interval{N}) where {N<:Real}
rand(::Type{Interval})
isflat(::Interval)
plot_recipe(::Interval{N}, ::N=zero(N)) where {N<:Real}
linear_map(::AbstractMatrix{N}, ::Interval{N}) where {N<:Real}
scale(::N, ::Interval{N}) where {N<:Real}
constraints_list(::Interval{N}) where {N<:Real}
center(::Interval{N}, ::Int) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N<:Real})
* [`singleton_list`](@ref singleton_list(::AbstractPolytope{N}) where {N<:Real})

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))

Inherited from [`AbstractZonotope`](@ref):
* [`ngens`](@ref ngens(::AbstractZonotope))
* [`order`](@ref order(::AbstractZonotope))
* [`togrep`](@ref togrep(::AbstractZonotope))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real))
* [`generators`](@ref generators(::AbstractHyperrectangle))
* [`genmat`](@ref genmat(::AbstractHyperrectangle))
