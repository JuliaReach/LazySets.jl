```@meta
CurrentModule = LazySets.HPolygonModule
```

# [Polygon in constraint representation (HPolygon)](@id def_HPolygon)

```@docs
HPolygon
```

## Operations

```@docs
σ(::AbstractVector, ::HPolygon)
translate(::HPolygon, ::AbstractVector)
```

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`reflect`](@ref reflect(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::LazySet)

Inherited from [`AbstractPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N})

Inherited from [`AbstractPolygon`](@ref):
* [`dim`](@ref dim(::AbstractPolygon))
* [`volume`](@ref volume(::AbstractPolygon))

Inherited from [`AbstractHPolygon`](@ref):
* [`an_element`](@ref an_element(::AbstractHPolygon))
* [`constraints_list`](@ref constraints_list(::AbstractHPolygon))
* [`isbounded`](@ref isbounded(::AbstractHPolygon, ::Bool=true))
* [`normalize`](@ref normalize(::AbstractHPolygon{N}, p::Real=N(2)) where {N})
* [`remove_redundant_constraints!`](@ref remove_redundant_constraints!(::AbstractHPolygon))
* [`tohrep`](@ref tohrep(::HPOLYGON) where {HPOLYGON<:AbstractHPolygon})
* [`tovrep`](@ref tovrep(::AbstractHPolygon))
* [`vertices_list`](@ref vertices_list(::AbstractHPolygon{N}) where {N})
* [`addconstraint!`](@ref addconstraint!(::AbstractHPolygon, ::HalfSpace))
* [`∈`](@ref ∈(::AbstractVector, ::AbstractHPolygon))
* [`isredundant`](@ref isredundant(::HalfSpace, ::HalfSpace, ::HalfSpace))
