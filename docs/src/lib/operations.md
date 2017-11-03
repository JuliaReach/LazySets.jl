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
dim(ms::MinkowskiSum)
MinkowskiSumArray
dim(msa::MinkowskiSumArray)
```

## Convex Hull

```@docs
ConvexHull
```

## Cartesian Product

```@docs
CartesianProduct
dim(cp::CartesianProduct)
CartesianProductArray
dim(cp::CartesianProductArray)
is_contained(d::Vector{Float64}, cp::CartesianProduct)
σ(d::AbstractVector{Float64}, cp::CartesianProductArray)
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