```@contents
Pages = ["intersection.md"]
Depth = 3
```

# Intersection of two sets

```@meta
CurrentModule = LazySets.API
```

```@docs; canonical=false
intersection(::LazySet, ::LazySet)
```

```@meta
CurrentModule = LazySets
```

```@docs
intersection(::AbstractHyperrectangle, ::AbstractHyperrectangle)
intersection(::Interval, ::HalfSpace)
intersection(::Interval, ::LazySet)
intersection(::AbstractHPolygon, ::AbstractHPolygon)
intersection(::Union{VPolygon,VPolytope}, ::Union{VPolygon,VPolytope})
intersection(::UnionSet, ::LazySet)
intersection(::AbstractPolyhedron, ::ResetMap)
intersection(::CartesianProductArray, ::CartesianProductArray)
intersection(::Union{CartesianProduct,CartesianProductArray}, ::AbstractPolyhedron)
intersection(::AbstractZonotope{N}, ::HalfSpace{N}) where {N}
intersection!(::Star, ::HalfSpace)
LazySets._bound_intersect_2D(::Zonotope, ::Line2D)
```
