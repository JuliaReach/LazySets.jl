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
* [`singleton_list`](@ref singleton_list(::LazySet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`linear_map`](@ref linear_map(::AbstractMatrix{NM}, ::AbstractPolyhedron{NP}) where {NM, NP})

Inherited from [`AbstractPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N})

Inherited from [`AbstractPolygon`](@ref):
* [`dim`](@ref dim(::AbstractPolygon))

Inherited from [`AbstractHPolygon`](@ref):
* [`an_element`](@ref an_element(::AbstractHPolygon))
* [`∈`](@ref ∈(::AbstractVector, ::AbstractHPolygon))
* [`vertices_list`](@ref vertices_list(::AbstractHPolygon{N}) where {N})
* [`tohrep`](@ref tohrep(::HPOLYGON) where {HPOLYGON<:AbstractHPolygon})
* [`tovrep`](@ref tovrep(::AbstractHPolygon))
* [`normalize`](@ref normalize(::AbstractHPolygon{N}, p=N(2)) where {N})
* [`isbounded`](@ref isbounded(::AbstractHPolygon, ::Bool=true))
* [`addconstraint!`](@ref addconstraint!(::AbstractHPolygon, ::LinearConstraint))
* [`isredundant`](@ref isredundant(::LinearConstraint, ::LinearConstraint, ::LinearConstraint))
* [`remove_redundant_constraints!`](@ref remove_redundant_constraints!(::AbstractHPolygon))
* [`constraints_list`](@ref constraints_list(::AbstractHPolygon))
