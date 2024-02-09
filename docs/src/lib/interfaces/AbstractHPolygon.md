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
∈(::AbstractVector, ::AbstractHPolygon)
rand(::Type{HPOLYGON}) where {HPOLYGON<:AbstractHPolygon}
tohrep(::HPOLYGON) where {HPOLYGON<:AbstractHPolygon}
tovrep(::AbstractHPolygon)
addconstraint!(::AbstractHPolygon, ::HalfSpace)
addconstraint!(::Vector{LC}, ::HalfSpace) where {LC<:HalfSpace}
normalize(P::AbstractHPolygon{N}, p::Real=N(2)) where {N}
isredundant(::HalfSpace, ::HalfSpace, ::HalfSpace)
remove_redundant_constraints!(::AbstractHPolygon)
constraints_list(::AbstractHPolygon)
vertices_list(::AbstractHPolygon{N}) where {N}
isbounded(::AbstractHPolygon, ::Bool=true)
```

## Implementations

* [Polygon in constraint representation (HPolygon)](@ref def_HPolygon)
* [Polygon in optimized constraint representation (HPolygonOpt)](@ref def_HPolygonOpt)
