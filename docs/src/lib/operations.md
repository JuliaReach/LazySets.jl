# Common Set Operations

This section of the manual describes the basic symbolic types describing
operations between sets.

```@contents
Pages = ["operations.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
DocTestSetup = quote
    using LazySets
end
```

## Cartesian Product

### Binary Cartesian Product

```@docs
CartesianProduct
×(::LazySet, ::LazySet)
*(::LazySet, ::LazySet)
dim(::CartesianProduct)
σ(::AbstractVector{Real}, ::CartesianProduct{Real})
∈(::AbstractVector{Real}, ::CartesianProduct{Real})
```

### ``n``-ary Cartesian Product

```@docs
CartesianProductArray
dim(::CartesianProductArray)
σ(::AbstractVector{Real}, ::CartesianProductArray{Real})
∈(::AbstractVector{Real}, ::CartesianProductArray{Real})
array(::CartesianProductArray)
```

## Convex Hull

### Binary Convex Hull

```@docs
ConvexHull
CH
dim(::ConvexHull)
σ(::AbstractVector{Real}, ::ConvexHull{Real})
```

### ``n``-ary Convex Hull

```@docs
ConvexHullArray
CHArray
dim(::ConvexHullArray)
σ(::AbstractVector{Real}, ::ConvexHullArray{Real})
array(::ConvexHullArray)
```

### Convex Hull Algorithms

```@docs
convex_hull
convex_hull!
right_turn
monotone_chain!
```

## Intersection

### Binary Intersection

```@docs
Intersection
∩(::LazySet, ::LazySet)
dim(::Intersection)
σ(::AbstractVector{Real}, ::Intersection{Real})
∈(::AbstractVector{Real}, ::Intersection{Real})
isempty(::Intersection)
```

### ``n``-ary Intersection

```@docs
IntersectionArray
dim(::IntersectionArray)
σ(::AbstractVector{Real}, ::IntersectionArray{Real})
array(::IntersectionArray)
```

## Minkowski Sum

### Binary Minkowski Sum

```@docs
MinkowskiSum
⊕(::LazySet, ::LazySet)
+(::LazySet, ::LazySet)
dim(::MinkowskiSum)
σ(::AbstractVector{Real}, ::MinkowskiSum{Real})
```

### ``n``-ary Minkowski Sum

```@docs
MinkowskiSumArray
dim(::MinkowskiSumArray)
σ(::AbstractVector{Real}, ::MinkowskiSumArray{Real})
array(::MinkowskiSumArray)
```

### ``n``-ary Minkowski Sum with cache

```@docs
CacheMinkowskiSum
dim(::CacheMinkowskiSum)
σ(::AbstractVector{Real}, ::CacheMinkowskiSum{Real})
array(::CacheMinkowskiSum)
forget_sets!(::CacheMinkowskiSum)
```

## Maps

### Linear Map

```@docs
LinearMap
*(::AbstractMatrix, ::LazySet)
*(::Real, ::LazySet)
dim(::LinearMap)
σ(::AbstractVector{Real}, ::LinearMap{Real, Real})
∈(::AbstractVector{Real}, ::LinearMap{Real, Real})
an_element(::LinearMap)
```

### Exponential Map

```@docs
ExponentialMap
dim(::ExponentialMap)
σ(::AbstractVector{Real}, ::ExponentialMap{Real})
∈(::AbstractVector{Real}, ::ExponentialMap{Real})
```

```@docs
ExponentialProjectionMap
dim(::ExponentialProjectionMap)
σ(::AbstractVector{Real}, ::ExponentialProjectionMap{Real})
```

```@docs
SparseMatrixExp
*(::SparseMatrixExp, ::LazySet)
```

```@docs
ProjectionSparseMatrixExp
*(::ProjectionSparseMatrixExp, ::LazySet)
```

## Symmetric Interval Hull

```@docs
SymmetricIntervalHull
dim(::SymmetricIntervalHull)
σ(::V, ::SymmetricIntervalHull{N}) where {N<:Real, V<:AbstractVector{N}}
an_element(::SymmetricIntervalHull{Float64, LazySet{Float64}})
```
