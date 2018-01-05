# Approximations

This section of the manual describes the Cartesian decomposition algorithms and
the approximation of high-dimensional convex sets using projections.

```@contents
Pages = ["approximations.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets.Approximations
DocTestSetup = quote
    using LazySets, LazySets.Approximations
end
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
LocalApproximation
PolygonalOverapproximation
new_approx(S::LazySet, p1::Vector{Float64}, d1::Vector{Float64}, p2::Vector{Float64}, d2::Vector{Float64})
addapproximation!(Ω::PolygonalOverapproximation, p1::Vector{Float64}, d1::Vector{Float64}, p2::Vector{Float64}, d2::Vector{Float64})
refine(Ω::PolygonalOverapproximation, i::Int)
tohrep(Ω::PolygonalOverapproximation)
approximate(S::LazySet, ɛ::Float64)
```

See [Iterative Refinement](@ref) for more details.
