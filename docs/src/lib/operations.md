# Common Set Operations

This section of the manual describes the basic symbolic types describing
operations between sets.

```@contents
Pages = ["operations.md"]
```

```@meta
CurrentModule = LazySets
```

## Minkowski Sum

```@docs
MinkowskiSum
MinkowskiSumArray
```

## Convex Hull

```@docs
ConvexHull
```

## Cartesian Product

```@docs
CartesianProduct
CartesianProductArray
```

## Linear Maps

```@docs
LinearMap
```

## Exponential Maps

```@docs
ExponentialMap
ExponentialProjectionMap
ProjectionSparseMatrixExp
SparseMatrixExp
σ(d::Union{Vector{Float64}, SparseVector{Float64,Int64}},
           eprojmap::ExponentialProjectionMap)
```