```@contents
Pages = ["minkowski_sum.md"]
Depth = 3
```

# Minkowski Sum

```@meta
CurrentModule = LazySets.API
```

```@docs; canonical=false
minkowski_sum(::LazySet, ::LazySet)
```

```@meta
CurrentModule = LazySets
```

```@docs
minkowski_sum(::LazySet, ::LazySet)
minkowski_sum(::AbstractPolyhedron, ::AbstractPolyhedron)
minkowski_sum(::AbstractHyperrectangle, ::AbstractHyperrectangle)
minkowski_sum(::AbstractZonotope, ::AbstractZonotope)
minkowski_sum(::DensePolynomialZonotope, ::AbstractZonotope)
minkowski_sum(::AbstractSingleton, ::AbstractSingleton)
```
