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

```@docs
Approximations
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
interval_hull
box_approximation_symmetric
symmetric_interval_hull
box_approximation_helper
```

### Metric properties of sets

```@docs
norm(::LazySet, ::Real=Inf)
radius(::LazySet, ::Real=Inf)
diameter(::LazySet, ::Real=Inf)
```

## Iterative refinement

```@docs
approximate
```

See [Iterative Refinement](@ref) for more details.
