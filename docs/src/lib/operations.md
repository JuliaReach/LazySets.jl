# Common Set Operations

This section of the manual describes the basic symbolic types describing
operations between sets.

```@contents
Pages = ["operations.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

## Minkowski Sum

```@docs
MinkowskiSum
dim(ms::MinkowskiSum)
σ(d::AbstractVector{Float64}, ms::MinkowskiSum)
```

```@docs
MinkowskiSumArray
dim(ms::MinkowskiSumArray)
σ(d::AbstractVector{Float64}, ms::MinkowskiSumArray)
Base.:+(msa::MinkowskiSumArray, sf::LazySet)
Base.:+(msa1::MinkowskiSumArray, msa2::MinkowskiSumArray)
Base.:+(msa::MinkowskiSumArray, sf::VoidSet)
```

## Convex Hull

```@docs
ConvexHull
dim(ch::ConvexHull)
σ(d::AbstractVector{Float64}, ch::ConvexHull)
```

### Convex Hull Algorithms

```@docs
convex_hull
convex_hull!
right_turn
monotone_chain!
```

## Cartesian Product

```@docs
CartesianProduct
dim(cp::CartesianProduct)
σ(d::AbstractVector{Float64}, cp::CartesianProduct)
is_contained(d::AbstractVector{Float64}, cp::CartesianProduct)
```

```@docs
CartesianProductArray
dim(cp::CartesianProductArray)
σ(d::AbstractVector{Float64}, cp::CartesianProductArray)
is_contained(d::AbstractVector{Float64}, cp::CartesianProductArray)
```

## Linear Maps

```@docs
LinearMap
dim(lm::LinearMap)
σ(d::AbstractVector{Float64}, lm::LinearMap)
```

## Exponential Maps

```@docs
ExponentialMap
dim(emap::ExponentialMap)
σ(d::AbstractVector{Float64}, eprojmap::ExponentialProjectionMap)
```

```@docs
ExponentialProjectionMap
dim(eprojmap::ExponentialProjectionMap)
```

```@docs
ProjectionSparseMatrixExp
SparseMatrixExp
```
