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

### Binary Minkowski Sum

```@docs
MinkowskiSum
dim(ms::MinkowskiSum)
σ(d::AbstractVector{Float64}, ms::MinkowskiSum)
Base.:+(::LazySet, ::LazySet)
```

### ``n``-ary Minkowski Sum

```@docs
MinkowskiSumArray
dim(ms::MinkowskiSumArray)
σ(d::AbstractVector{Float64}, ms::MinkowskiSumArray)
Base.:+(::MinkowskiSumArray, ::LazySet)
Base.:+(::LazySet, ::MinkowskiSumArray)
Base.:+(::MinkowskiSumArray, ::MinkowskiSumArray)
Base.:+(::MinkowskiSumArray, ::VoidSet)
```

## Cartesian Product

### Binary Cartesian Product

```@docs
CartesianProduct
dim(cp::CartesianProduct)
σ(d::AbstractVector{Float64}, cp::CartesianProduct)
Base.:*(::LazySet, ::LazySet)
is_contained(d::AbstractVector{Float64}, cp::CartesianProduct)
```

### ``n``-ary Cartesian Product

```@docs
CartesianProductArray
dim(cp::CartesianProductArray)
σ(d::AbstractVector{Float64}, cp::CartesianProductArray)
Base.:*(::CartesianProductArray, ::LazySet)
Base.:*(::LazySet, ::CartesianProductArray)
Base.:*(::CartesianProductArray, ::CartesianProductArray)
is_contained(d::AbstractVector{Float64}, cp::CartesianProductArray)
```

## Maps

### Linear Map

```@docs
LinearMap
dim(lm::LinearMap)
σ(d::AbstractVector{Float64}, lm::LinearMap)
```

### Exponential Map

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
