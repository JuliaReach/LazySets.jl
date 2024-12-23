```@contents
Pages = ["AbstractHPolygon.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# [Polygons in constraint representation (AbstractHPolygon)](@id def_AbstractHPolygon)

An HPolygon is a polygon in H-representation (or constraint representation).

```@docs
AbstractHPolygon
```

This interface defines the following functions:

```@docs
an_element(::AbstractHPolygon)
constraints_list(::AbstractHPolygon)
isbounded(::AbstractHPolygon, ::Bool=true)
normalize(P::AbstractHPolygon{N}, p::Real=N(2)) where {N}
rand(::Type{HPOLYGON}) where {HPOLYGON<:AbstractHPolygon}
remove_redundant_constraints!(::AbstractHPolygon)
tohrep(::HPOLYGON) where {HPOLYGON<:AbstractHPolygon}
tovrep(::AbstractHPolygon)
vertices_list(::AbstractHPolygon)
addconstraint!(::AbstractHPolygon, ::HalfSpace)
addconstraint!(::Vector{LC}, ::HalfSpace) where {LC<:HalfSpace}
∈(::AbstractVector, ::AbstractHPolygon)
isredundant(::HalfSpace, ::HalfSpace, ::HalfSpace)
```

## Implementations

* [Polygon in constraint representation (HPolygon)](@ref def_HPolygon)
* [Polygon in optimized constraint representation (HPolygonOpt)](@ref def_HPolygonOpt)
