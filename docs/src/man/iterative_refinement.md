# Iterative Refinement

This section of the manual describes the Cartesian decomposition algorithms and
the approximation of high-dimensional convex sets using projections.

[under construction]

```@contents
Pages = ["iterative_refinement.md"]
```

```@meta
CurrentModule = LazySets.Approximations
```

```@docs
Approximation2D
Approximation2D(::Vector{Float64}, ::Vector{Float64}, ::Vector{Float64}, ::Vector{Float64})
```

```@docs
refine(::LazySet, ::Approximation2D)
```
