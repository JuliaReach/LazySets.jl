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
```

### Convenience functions

```@docs
uniform_partition
```

## Overapproximations

```@docs
overapproximate
LazySets.Approximations._overapproximate_zonotope_vrep
LazySets.Approximations._overapproximate_zonotope_cpa
```

## Underapproximations

```@docs
underapproximate
```

## Approximations

```@docs
approximate
```

## Box Approximations

```@docs
ballinf_approximation
box_approximation
interval_hull
symmetric_interval_hull
box_approximation_symmetric
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
_approximate(S::LazySet{N}, ε::Real) where {N<:AbstractFloat}
constraint(::LocalApproximation)
```

## Template directions

```@docs
AbstractDirections
isbounding
isnormalized
project(::LazySet, ::AbstractVector{Int}, ::Type{<:AbstractDirections})
BoxDirections
DiagDirections
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
