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
dim(::MinkowskiSum)
σ(::AbstractVector{Float64}, ::MinkowskiSum)
Base.:+(::LazySet, ::LazySet)
```

### ``n``-ary Minkowski Sum

```@docs
MinkowskiSumArray
dim(::MinkowskiSumArray)
σ(::AbstractVector{Float64}, ::MinkowskiSumArray)
Base.:+(::MinkowskiSumArray, ::LazySet)
Base.:+(::LazySet, ::MinkowskiSumArray)
Base.:+(::MinkowskiSumArray, ::MinkowskiSumArray)
Base.:+(::MinkowskiSumArray, ::VoidSet)
```

## Cartesian Product

### Binary Cartesian Product

```@docs
CartesianProduct
dim(::CartesianProduct)
σ(::AbstractVector{Float64}, ::CartesianProduct)
Base.:*(::LazySet, ::LazySet)
is_contained(::AbstractVector{Float64}, ::CartesianProduct)
```

### ``n``-ary Cartesian Product

```@docs
CartesianProductArray
dim(::CartesianProductArray)
σ(::AbstractVector{Float64}, ::CartesianProductArray)
Base.:*(::CartesianProductArray, ::LazySet)
Base.:*(::LazySet, ::CartesianProductArray)
Base.:*(::CartesianProductArray, ::CartesianProductArray)
is_contained(::AbstractVector{Float64}, ::CartesianProductArray)
```

## Maps

### Linear Map

```@docs
LinearMap
dim(::LinearMap)
σ(::AbstractVector{Float64}, ::LinearMap)
*(::AbstractMatrix{Float64}, ::LazySet)
*(::Real, ::LazySet)
```

### Exponential Map

```@docs
ExponentialMap
dim(::ExponentialMap)
σ(::AbstractVector{Float64}, ::ExponentialMap)
```

```@docs
ExponentialProjectionMap
dim(::ExponentialProjectionMap)
σ(::AbstractVector{Float64}, ::ExponentialProjectionMap)
```

```@docs
SparseMatrixExp
*(::SparseMatrixExp, ::LazySet)
```

```@docs
ProjectionSparseMatrixExp
*(::ProjectionSparseMatrixExp, ::LazySet)
```

## Convex Hull

```@docs
ConvexHull
dim(::ConvexHull)
σ(::AbstractVector{Float64}, ::ConvexHull)
```

### Convex Hull Algorithms

```@docs
convex_hull
convex_hull!
right_turn
monotone_chain!
```
