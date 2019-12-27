```@meta
CurrentModule = LazySets
```

# [Polygon in constraint representation (HPolygon)](@id def_HPolygon)

```@docs
HPolygon
σ(::AbstractVector{N}, ::HPolygon{N}) where {N<:Real}
translate(::HPolygon{N}, ::AbstractVector{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N<:Real})
* [`singleton_list`](@ref singleton_list(::AbstractPolytope{N}) where {N<:Real})
* [`linear_map`](@ref linear_map(::AbstractMatrix{N}, ::AbstractPolyhedron{N}) where {N<:Real})

Inherited from [`AbstractPolygon`](@ref):
* [`dim`](@ref dim(::AbstractPolygon))

Inherited from [`AbstractHPolygon`](@ref):
* [`an_element`](@ref an_element(::AbstractHPolygon{N}) where {N<:Real})
* [`∈`](@ref ∈(::AbstractVector{N}, ::AbstractHPolygon{N}) where {N<:Real})
* [`vertices_list`](@ref vertices_list(::AbstractHPolygon{N}) where {N<:Real})
* [`tohrep`](@ref tohrep(::HPOLYGON) where {HPOLYGON<:AbstractHPolygon})
* [`tovrep`](@ref tovrep(::AbstractHPolygon{N}) where {N<:Real})
* [`isbounded`](@ref isbounded(::AbstractHPolygon, ::Bool=true))
* [`addconstraint!`](@ref addconstraint!(::AbstractHPolygon{N}, ::LinearConstraint{N}) where {N<:Real})
* [`isredundant`](@ref isredundant(::LinearConstraint{N}, ::LinearConstraint{N}, ::LinearConstraint{N}) where {N<:Real})
* [`remove_redundant_constraints!`](@ref remove_redundant_constraints!(::AbstractHPolygon))
* [`constraints_list`](@ref constraints_list(::AbstractHPolygon{N}) where {N<:Real})
