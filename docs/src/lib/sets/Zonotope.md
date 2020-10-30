```@meta
CurrentModule = LazySets
```

# [Zonotope](@id def_Zonotope)

```@docs
Zonotope
center(::Zonotope{N}) where {N<:Real}
rand(::Type{Zonotope})
generators(::Zonotope)
genmat(::Zonotope)
scale(::Real, ::Zonotope)
scale!(::Real, Z::Zonotope)
ngens(::Zonotope)
togrep(::Zonotope)
reduce_order(::Zonotope, ::Union{Integer, Rational})
split(::AbstractZonotope, ::Int)
split(::AbstractZonotope, ::AbstractVector{Int}, ::AbstractVector{Int})
remove_zero_generators(::Zonotope)
linear_map!(::Zonotope, ::AbstractMatrix, ::Zonotope)
quadratic_map(::Vector{MT}, ::Zonotope{N}) where {N, MT<:AbstractMatrix{N}}
LazySets._bound_intersect_2D(::Zonotope, ::Line2D)
```

Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`singleton_list`](@ref singleton_list(::LazySet))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N<:Real})

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope{N}) where {N<:Real})

Inherited from [`AbstractZonotope`](@ref):
* [`ρ`](@ref ρ(::AbstractVector{N}, ::AbstractZonotope{N}) where {N<:Real})
* [`σ`](@ref σ(::AbstractVector{N}, ::AbstractZonotope{N}) where {N<:Real})
* [`∈`](@ref ∈(::AbstractVector{N}, ::AbstractZonotope{N}) where {N<:Real})
* [`linear_map`](@ref linear_map(::AbstractMatrix{N}, ::AbstractZonotope{N}) where {N<:Real})
* [`translate`](@ref translate(::AbstractZonotope{N}, ::AbstractVector{N}) where {N<:Real})
* [`constraints_list`](@ref constraints_list(::AbstractZonotope{N}) where {N<:Real})
* [`constraints_list`](@ref constraints_list(::AbstractZonotope{N}; ::Bool=true) where {N<:AbstractFloat})
* [`vertices_list`](@ref vertices_list(::AbstractZonotope{N}) where {N<:Real})
* [`order`](@ref order(::AbstractZonotope))
