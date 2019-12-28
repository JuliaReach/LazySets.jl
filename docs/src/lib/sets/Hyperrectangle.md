```@meta
CurrentModule = LazySets
```

# [Hyperrectangle](@id def_Hyperrectangle)

```@docs
Hyperrectangle
rand(::Type{Hyperrectangle})
center(::Hyperrectangle{N}) where {N<:Real}
radius_hyperrectangle(::Hyperrectangle{N}) where {N<:Real}
radius_hyperrectangle(::Hyperrectangle{N}, ::Int) where {N<:Real}
translate(::Hyperrectangle{N}, ::AbstractVector{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N<:Real})
* [`singleton_list`](@ref singleton_list(::AbstractPolytope{N}) where {N<:Real})

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope{N}) where {N<:Real})

Inherited from [`AbstractZonotope`](@ref):
* [`ngens`](@ref ngens(::AbstractZonotope))
* [`linear_map`](@ref linear_map(::AbstractMatrix{N}, ::AbstractZonotope{N}) where {N<:Real})
* [`order`](@ref order(::AbstractZonotope))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`σ`](@ref σ(::AbstractVector{N}, ::AbstractHyperrectangle{N}) where {N<:Real})
* [`ρ`](@ref ρ(::AbstractVector{N}, ::AbstractHyperrectangle{N}) where {N<:Real})
* [`∈`](@ref ∈(::AbstractVector{N}, ::AbstractHyperrectangle{N}) where {N<:Real})
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real))
* [`vertices_list`](@ref vertices_list(::AbstractHyperrectangle{N}) where {N<:Real})
* [`high`](@ref high(::AbstractHyperrectangle{N}) where {N<:Real})
* [`low`](@ref low(::AbstractHyperrectangle{N}) where {N<:Real})
* [`isflat`](@ref isflat(::Hyperrectangle))
* [`generators`](@ref generators(::AbstractHyperrectangle))
* [`genmat`](@ref genmat(::AbstractHyperrectangle{N}) where {N<:Real})
* [`constraints_list`](@ref constraints_list(::AbstractHyperrectangle{N}) where {N<:Real})
