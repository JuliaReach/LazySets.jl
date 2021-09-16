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
low(::Interval{N}) where {N}
high(::Interval{N}) where {N}
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
```
Inherited from [`LazySet`](@ref):
* [`singleton_list`](@ref singleton_list(::LazySet))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N})

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
