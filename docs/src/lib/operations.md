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
    using Compat.SparseArrays, Compat.LinearAlgebra
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
vertices_list(::CartesianProduct{Real})
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{Real}))

### ``n``-ary Cartesian Product

```@docs
CartesianProductArray
dim(::CartesianProductArray)
σ(::AbstractVector{Real}, ::CartesianProductArray{Real})
∈(::AbstractVector{Real}, ::CartesianProductArray{Real})
array(::CartesianProductArray)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{Real}))

## Convex Hull

### Binary Convex Hull

```@docs
ConvexHull
CH
dim(::ConvexHull)
σ(::AbstractVector{Real}, ::ConvexHull{Real})
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{Real}))

### ``n``-ary Convex Hull

```@docs
ConvexHullArray
CHArray
dim(::ConvexHullArray)
σ(::AbstractVector{Real}, ::ConvexHullArray{Real})
array(::ConvexHullArray)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{Real}))

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
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{Real}))

### ``n``-ary Intersection

```@docs
IntersectionArray
dim(::IntersectionArray)
σ(::AbstractVector{Real}, ::IntersectionArray{Real})
array(::IntersectionArray)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{Real}))

## Minkowski Sum

### Binary Minkowski Sum

```@docs
MinkowskiSum
⊕(::LazySet, ::LazySet)
+(::LazySet, ::LazySet)
dim(::MinkowskiSum)
σ(::AbstractVector{Real}, ::MinkowskiSum{Real})
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{Real}))

### ``n``-ary Minkowski Sum

```@docs
MinkowskiSumArray
dim(::MinkowskiSumArray)
σ(::AbstractVector{Real}, ::MinkowskiSumArray{Real})
array(::MinkowskiSumArray)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{Real}))

### ``n``-ary Minkowski Sum with cache

```@docs
CacheMinkowskiSum
dim(::CacheMinkowskiSum)
σ(::AbstractVector{Real}, ::CacheMinkowskiSum{Real})
array(::CacheMinkowskiSum)
forget_sets!(::CacheMinkowskiSum)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{Real}))

## Maps

### Linear Map

```@docs
LinearMap
*(::AbstractMatrix{Real}, ::LazySet{Real})
*(::Real, ::LazySet{Real})
dim(::LinearMap)
σ(::AbstractVector{Real}, ::LinearMap{Real})
∈(::AbstractVector{Real}, ::LinearMap{Real, LazySet{Real}, Real, Matrix{Real}})
an_element(::LinearMap)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

### Exponential Map

```@docs
ExponentialMap
dim(::ExponentialMap)
σ(::AbstractVector{Real}, ::ExponentialMap{Real})
∈(::AbstractVector{Real}, ::ExponentialMap{Real})
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{Real}))

```@docs
ExponentialProjectionMap
dim(::ExponentialProjectionMap)
σ(::AbstractVector{Real}, ::ExponentialProjectionMap{Real})
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{Real}))

```@docs
SparseMatrixExp
*(::SparseMatrixExp, ::LazySet)
get_row(::SparseMatrixExp, ::Int)
```

```@docs
ProjectionSparseMatrixExp
*(::ProjectionSparseMatrixExp, ::LazySet)
```

## Symmetric Interval Hull

```@docs
SymmetricIntervalHull
σ(::AbstractVector{N}, ::SymmetricIntervalHull{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`singleton_list`](@ref singleton_list(::AbstractPolytope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractPolytope))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* `an_element`

Inherited from [`AbstractHyperrectangle`](@ref):
* [`∈`](@ref ∈(::AbstractVector{Real}, ::AbstractHyperrectangle{Real}))
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real))
* [`vertices_list`](@ref vertices_list(::AbstractHyperrectangle{Real}))
* [`high`](@ref high(::AbstractHyperrectangle{Real}))
* [`low`](@ref low(::AbstractHyperrectangle{Real}))
