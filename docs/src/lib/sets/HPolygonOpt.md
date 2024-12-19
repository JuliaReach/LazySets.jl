```@meta
CurrentModule = LazySets
```

# [Polygon in optimized constraint representation (HPolygonOpt)](@id def_HPolygonOpt)

```@docs
HPolygonOpt
σ(::AbstractVector, ::HPolygonOpt)
translate(::HPolygonOpt, ::AbstractVector)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`reflect`](@ref reflect(::LazySet))

Inherited from [`AbstractPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope, ::Bool=false))

Inherited from [`AbstractPolygon`](@ref):
* [`dim`](@ref dim(::AbstractPolygon))
* [`volume`](@ref volume(::AbstractPolygon))

Inherited from [`AbstractHPolygon`](@ref):
* [`an_element`](@ref an_element(::AbstractHPolygon))
* [`∈`](@ref ∈(::AbstractVector, ::AbstractHPolygon))
* [`vertices_list`](@ref vertices_list(::AbstractHPolygon{N}) where {N})
* [`tohrep`](@ref tohrep(::HPOLYGON) where {HPOLYGON<:AbstractHPolygon})
* [`tovrep`](@ref tovrep(::AbstractHPolygon))
* [`normalize`](@ref normalize(::AbstractHPolygon{N}, p::Real=N(2)) where {N})
* [`isbounded`](@ref isbounded(::AbstractHPolygon, ::Bool=true))
* [`addconstraint!`](@ref addconstraint!(::AbstractHPolygon, ::HalfSpace))
* [`isredundant`](@ref isredundant(::HalfSpace, ::HalfSpace, ::HalfSpace))
* [`remove_redundant_constraints!`](@ref remove_redundant_constraints!(::AbstractHPolygon))
* [`constraints_list`](@ref constraints_list(::AbstractHPolygon))
