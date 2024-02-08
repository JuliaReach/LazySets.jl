```@contents
Pages = ["issubset.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# Subset Check

```@docs
issubset
⊆(::LazySet, ::LazySet, ::Bool=false)
⊆(::LazySet, ::AbstractHyperrectangle, ::Bool=false)
⊆(::AbstractPolytope, ::LazySet, ::Bool=false)
⊆(::AbstractPolytope, ::AbstractHyperrectangle, ::Bool=false)
⊆(::AbstractZonotope, ::AbstractHyperrectangle)
⊆(::AbstractHyperrectangle, ::AbstractHyperrectangle, ::Bool=false)
⊆(::LazySet, ::AbstractPolyhedron, ::Bool=false)
⊆(::AbstractSingleton, ::LazySet, ::Bool=false)
⊆(::AbstractSingleton, ::AbstractHyperrectangle, ::Bool=false)
⊆(::AbstractSingleton, ::AbstractSingleton, ::Bool=false)
⊆(::Ball2, ::Ball2, ::Bool=false)
⊆(::Union{Ball2, Ballp}, ::AbstractSingleton, ::Bool=false)
⊆(::LineSegment, ::LazySet, ::Bool=false)
⊆(::LineSegment, ::AbstractHyperrectangle, ::Bool=false)
⊆(::Interval, ::Interval, ::Bool=false)
⊆(::Interval, ::UnionSet, ::Bool=false)
⊆(::EmptySet, ::LazySet, ::Bool=false)
⊆(::LazySet, ::EmptySet, ::Bool=false)
⊆(::UnionSet, ::LazySet, ::Bool=false)
⊆(::UnionSetArray, ::LazySet, ::Bool=false)
⊆(::LazySet, ::Universe, ::Bool=false)
⊆(::Universe, ::LazySet, ::Bool=false)
⊆(::LazySet, ::Complement, ::Bool=false)
⊆(::CartesianProduct, ::CartesianProduct, ::Bool=false)
⊆(::CartesianProductArray, ::CartesianProductArray, ::Bool=false)
⊆(::AbstractZonotope, ::AbstractHyperrectangle, ::Bool=false)
⊆(::LazySet, ::UnionSetArray, ::Bool=false; ::Bool=true)
⊂
```
