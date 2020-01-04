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
project
```

### Convenience functions

```@docs
uniform_partition
```

## Overapproximations

```@docs
overapproximate
```

## Underapproximations

```@docs
underapproximate
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

## Iterative refinement

```@docs
LocalApproximation
PolygonalOverapproximation
new_approx(S::LazySet, p1::VN, d1::VN,
           p2::VN, d2::VN) where {N<:AbstractFloat, VN<:AbstractVector{N}}
addapproximation!(Ω::PolygonalOverapproximation, p1::VN, d1::VN, p2::VN, d2::VN) where {N<:Real, VN<:AbstractVector{N}}
refine(::LocalApproximation, ::LazySet)
tohrep(::PolygonalOverapproximation)
approximate(S::LazySet{N}, ε::N) where {N<:AbstractFloat}
constraint(::LocalApproximation)
```

## Template directions

```@docs
AbstractDirections
isbounding
BoxDirections
OctDirections
BoxDiagDirections
PolarDirections
SphericalDirections
CustomDirections
```

See also `overapproximate(X::LazySet, dir::AbstractDirections)::HPolytope`.

## Hausdorff distance

```@docs
hausdorff_distance
```
