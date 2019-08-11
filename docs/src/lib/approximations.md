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
new_approx(::LazySet, ::Vector{Float64}, ::Vector{Float64}, ::Vector{Float64}, ::Vector{Float64})
addapproximation!(::PolygonalOverapproximation, ::Vector{Float64}, ::Vector{Float64}, ::Vector{Float64}, ::Vector{Float64})
refine(::LocalApproximation, ::LazySet)
tohrep(::PolygonalOverapproximation)
approximate(::LazySet{Float64}, ::Float64)
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
