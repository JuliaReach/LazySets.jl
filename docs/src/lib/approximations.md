# Approximations

This section of the manual describes the Cartesian decomposition algorithms and
the approximation of high-dimensional convex sets using projections.

```@contents
Pages = ["approximations.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets.Approximations
```

## Cartesian Decomposition

```@docs
decompose
overapproximate
```

## Box Approximations

```@docs
ballinf_approximation
box_approximation
box_approximation_symmetric
box_approximation_helper
```

### Metric properties of sets

```@docs
norm(X::LazySet, p=Inf)
radius(X::LazySet, p=Inf)
diameter(X::LazySet, p=Inf)
```

## Iterative refinement


```@docs
approximate
```
