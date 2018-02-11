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
Base.:*(::LazySet{Float64}, ::LazySet{Float64})
dim(::CartesianProduct{Float64, LazySet{Float64}, LazySet{Float64}})
σ(::AbstractVector{Float64}, ::CartesianProduct{Float64, LazySet{Float64}, LazySet{Float64}})
∈(::AbstractVector{Float64}, ::CartesianProduct{Float64, LazySet{Float64}, LazySet{Float64}})
```

### ``n``-ary Cartesian Product

```@docs
CartesianProductArray{Float64, LazySet{Float64}}
array(::CartesianProductArray{Float64, LazySet{Float64}})
dim(::CartesianProductArray{Float64, LazySet{Float64}})
σ(::AbstractVector{Float64}, ::CartesianProductArray{Float64, LazySet{Float64}})
∈(::AbstractVector{Float64}, ::CartesianProductArray{Float64, LazySet{Float64}})
```

## Convex Hull

### Binary Convex Hull

```@docs
ConvexHull
CH
dim(::ConvexHull{Float64, LazySet{Float64}, LazySet{Float64}})
σ(::AbstractVector{Float64}, ::ConvexHull{Float64, LazySet{Float64}, LazySet{Float64}})
```

### ``n``-ary Convex Hull

```@docs
ConvexHullArray
CHArray
array(::ConvexHullArray{Float64, LazySet{Float64}})
dim(cha::ConvexHullArray)
σ(d::AbstractVector{Float64}, cha::ConvexHullArray)
```

### Convex Hull Algorithms

```@docs
convex_hull
convex_hull!
right_turn
monotone_chain!
```

## Intersection

```@docs
Intersection
dim(::Intersection{Float64, LazySet{Float64}, LazySet{Float64}})
σ(::AbstractVector{Float64}, ::Intersection{Float64, LazySet{Float64}, LazySet{Float64}})
∈(::AbstractVector{Float64}, ::Intersection{Float64, LazySet{Float64}, LazySet{Float64}})
isempty(::Intersection{Float64, LazySet{Float64}, LazySet{Float64}})
```

## Minkowski Sum

### Binary Minkowski Sum

```@docs
MinkowskiSum
Base.:+(::LazySet{Float64}, ::LazySet{Float64})
⊕
dim(::MinkowskiSum{Float64, LazySet{Float64}, LazySet{Float64}})
σ(::AbstractVector{Float64}, ::MinkowskiSum{Float64, LazySet{Float64}, LazySet{Float64}})
```

### ``n``-ary Minkowski Sum

```@docs
MinkowskiSumArray
array(::MinkowskiSumArray{Float64, LazySet{Float64}})
dim(::MinkowskiSumArray{Float64, LazySet{Float64}})
σ(::AbstractVector{Float64}, ::MinkowskiSumArray{Float64, LazySet{Float64}})
```

## Maps

### Linear Map

```@docs
LinearMap
dim(::LinearMap{Float64, Float64})
σ(::AbstractVector{Float64}, ::LinearMap{Float64, Float64})
*(::AbstractMatrix, ::LazySet)
*(::Float64, ::LazySet)
∈(x::AbstractVector{Float64}, ::LinearMap{Float64, Float64})
```

### Exponential Map

```@docs
ExponentialMap
dim(::ExponentialMap{Float64, LazySet{Float64}})
σ(::AbstractVector{Float64}, ::ExponentialMap{Float64, LazySet{Float64}})
∈(::AbstractVector{Float64}, ::ExponentialMap{Float64, LazySet{Float64}})
```

```@docs
ExponentialProjectionMap
dim(::ExponentialProjectionMap{Float64, LazySet{Float64}})
σ(::AbstractVector{Float64}, ::ExponentialProjectionMap{Float64, LazySet{Float64}})
```

```@docs
SparseMatrixExp
*(::SparseMatrixExp{Float64}, ::LazySet{Float64})
```

```@docs
ProjectionSparseMatrixExp
*(::ProjectionSparseMatrixExp{Float64}, ::LazySet{Float64})
```

## Symmetric Interval Hull

```@docs
SymmetricIntervalHull
dim(::SymmetricIntervalHull{Float64, LazySet{Float64}})
σ(::AbstractVector{Float64}, ::SymmetricIntervalHull{Float64, LazySet{Float64}})
an_element(::SymmetricIntervalHull{Float64, LazySet{Float64}})
```
