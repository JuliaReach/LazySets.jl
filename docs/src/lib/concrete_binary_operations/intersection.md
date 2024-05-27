```@contents
Pages = ["intersection.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# Intersection of two sets

```@docs
intersection(::AbstractSingleton, ::LazySet)
intersection(::Line2D, ::Line2D)
intersection(::AbstractHyperrectangle, ::AbstractHyperrectangle)
intersection(::Interval, ::HalfSpace)
intersection(::Interval, ::Hyperplane)
intersection(::Interval, ::LazySet)
intersection(::AbstractHPolygon, ::AbstractHPolygon)
intersection(::AbstractPolyhedron{N}, ::AbstractPolyhedron{N}) where {N}
intersection(::Union{VPolytope, VPolygon}, ::Union{VPolytope, VPolygon})
intersection(::VPolygon, ::VPolygon; ::Bool=true)
intersection(::UnionSet, ::LazySet)
intersection(::UnionSetArray, ::LazySet)
intersection(::Universe, ::LazySet)
intersection(::AbstractPolyhedron, ::ResetMap)
intersection(::CartesianProductArray, ::CartesianProductArray)
intersection(::LinearMap, ::LazySet)
intersection(::Union{CartesianProduct,CartesianProductArray}, ::AbstractPolyhedron)
intersection(::LineSegment, ::Line2D)
intersection(::LineSegment, ::LineSegment)
intersection(::AbstractZonotope{N}, ::HalfSpace{N}) where {N}
intersection(::Star, ::HalfSpace)
intersection!(::Star, ::HalfSpace)
```
