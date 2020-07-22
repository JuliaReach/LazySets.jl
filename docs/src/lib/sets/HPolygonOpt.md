```@meta
CurrentModule = LazySets
```

# [Polygon in optimized constraint representation (HPolygonOpt)](@id def_HPolygonOpt)

```@docs
HPolygonOpt
σ(::AbstractVector{N}, ::HPolygonOpt{N}) where {N<:Real}
translate(::HPolygonOpt{N}, ::AbstractVector{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`singleton_list`](@ref singleton_list(::LazySet))

Inherited from [`AbstractPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N<:Real})

Inherited from [`AbstractPolygon`](@ref):
* [`dim`](@ref dim(::AbstractPolygon))

Inherited from [`AbstractHPolygon`](@ref):
* [`an_element`](@ref an_element(::AbstractHPolygon{N}) where {N<:Real})
* [`∈`](@ref ∈(::AbstractVector{N}, ::AbstractHPolygon{N}) where {N<:Real})
* [`vertices_list`](@ref vertices_list(::AbstractHPolygon{N}) where {N<:Real})
* [`tohrep`](@ref tohrep(::HPOLYGON) where {HPOLYGON<:AbstractHPolygon})
* [`tovrep`](@ref tovrep(::AbstractHPolygon{N}) where {N<:Real})
* [`normalize`](@ref normalize(::AbstractHPolygon{N}, p=N(2)) where {N<:Real})
* [`isbounded`](@ref isbounded(::AbstractHPolygon, ::Bool=true))
* [`addconstraint!`](@ref addconstraint!(::AbstractHPolygon{N}, ::LinearConstraint{N}) where {N<:Real})
* [`isredundant`](@ref isredundant(::LinearConstraint{N}, ::LinearConstraint{N}, ::LinearConstraint{N}) where {N<:Real})
* [`remove_redundant_constraints!`](@ref remove_redundant_constraints!(::AbstractHPolygon))
* [`constraints_list`](@ref constraints_list(::AbstractHPolygon{N}) where {N<:Real})
