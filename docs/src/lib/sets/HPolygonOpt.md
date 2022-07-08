```@meta
CurrentModule = LazySets
```

# [Polygon in optimized constraint representation (HPolygonOpt)](@id def_HPolygonOpt)

```@docs
HPolygonOpt
σ(::AbstractVector, ::HPolygonOpt)
translate(::HPolygonOpt, ::AbstractVector)
```
Inherited from [`ConvexSet`](@ref):
* [`norm`](@ref norm(::ConvexSet, ::Real))
* [`radius`](@ref radius(::ConvexSet, ::Real))
* [`diameter`](@ref diameter(::ConvexSet, ::Real))
* [`singleton_list`](@ref singleton_list(::ConvexSet))

Inherited from [`AbstractPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N})
* [`volume`](@ref volume(::AbstractPolytope))

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
