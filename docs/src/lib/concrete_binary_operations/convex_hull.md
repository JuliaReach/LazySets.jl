```@contents
Pages = ["convex_hull.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# Convex Hull

```@docs
convex_hull(::LazySet, ::LazySet)
convex_hull(::HPoly, ::HPoly)
convex_hull(::Vector{VN}) where {N, VN<:AbstractVector{N}}
monotone_chain!
```
