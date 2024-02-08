```@contents
Pages = ["isdisjoint.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# Check for Disjointness of Sets

!!! note
    `is_intersection_empty` can be used as an alternative name to `isdisjoint`.

```@docs
isdisjoint(::LazySet, ::LazySet, ::Bool=false)
isdisjoint(::AbstractHyperrectangle, ::AbstractHyperrectangle, ::Bool=false)
isdisjoint(::LazySet, ::AbstractSingleton, ::Bool=false)
isdisjoint(::AbstractSingleton, ::AbstractSingleton, ::Bool=false)
isdisjoint(::AbstractZonotope, ::Hyperplane, ::Bool=false)
isdisjoint(::Ball2, ::Ball2, ::Bool=false)
isdisjoint(::LineSegment, ::LineSegment, ::Bool=false)
isdisjoint(::LazySet, ::Hyperplane, ::Bool=false)
isdisjoint(::LazySet, ::HalfSpace, ::Bool=false)
isdisjoint(::HalfSpace, ::HalfSpace, ::Bool=false)
isdisjoint(::AbstractPolyhedron, ::LazySet, ::Bool=false)
isdisjoint(::UnionSet, ::LazySet, ::Bool=false)
isdisjoint(::UnionSetArray, ::LazySet, ::Bool=false)
isdisjoint(::Universe, ::LazySet, ::Bool=false)
isdisjoint(::Complement, ::LazySet, ::Bool=false)
isdisjoint(::AbstractZonotope, ::AbstractZonotope, ::Bool=false)
isdisjoint(::Interval, ::Interval, ::Bool=false)
isdisjoint(::CartesianProductArray, ::AbstractPolyhedron, ::Bool=false)
isdisjoint(::CartesianProductArray, ::CartesianProductArray, ::Bool=false)
isdisjoint(::CartesianProductArray, ::AbstractHyperrectangle, ::Bool=false)
isdisjoint(::Line2D, ::Line2D, ::Bool=false)
```
