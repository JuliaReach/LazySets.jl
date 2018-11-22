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
ρ(::AbstractVector{N}, ::CartesianProduct{N}) where {N<:Real}
σ(::AbstractVector{N}, ::CartesianProduct{N}) where {N<:Real}
∈(::AbstractVector{N}, ::CartesianProduct{N}) where {N<:Real}
isempty(::CartesianProduct)
constraints_list(::CartesianProduct{N}) where {N<:Real}
vertices_list(::CartesianProduct{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{N}) where {N<:Real})

### ``n``-ary Cartesian Product

```@docs
CartesianProductArray
dim(::CartesianProductArray)
ρ(::AbstractVector{N}, ::CartesianProductArray{N}) where {N<:Real}
σ(::AbstractVector{N}, ::CartesianProductArray{N}) where {N<:Real}
∈(::AbstractVector{N}, ::CartesianProductArray{N}) where {N<:Real}
isempty(::CartesianProductArray)
constraints_list(::CartesianProductArray{N}) where {N<:Real}
vertices_list(::CartesianProductArray{N}) where {N<:Real}
array(::CartesianProductArray{N, S}) where {N<:Real, S<:LazySet{N}}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{N}) where {N<:Real})

## Convex Hull

### Binary Convex Hull

```@docs
ConvexHull
CH
dim(::ConvexHull)
ρ(::AbstractVector{N}, ::ConvexHull{N}) where {N<:Real}
σ(::AbstractVector{N}, ::ConvexHull{N}) where {N<:Real}
isempty(::ConvexHull)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{N}) where {N<:Real})

### ``n``-ary Convex Hull

```@docs
ConvexHullArray
CHArray
dim(::ConvexHullArray)
ρ(::AbstractVector{N}, ::ConvexHullArray{N}) where {N<:Real}
σ(::AbstractVector{N}, ::ConvexHullArray{N}) where {N<:Real}
array(::ConvexHullArray{N, S}) where {N<:Real, S<:LazySet{N}}
isempty(::ConvexHullArray)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{N}) where {N<:Real})

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
ρ(::AbstractVector{N}, ::Intersection{N}) where {N<:Real}
ρ(::AbstractVector{N}, ::Intersection{N, S1, S2}) where {N<:Real, S1<:LazySet{N}, S2<:Union{HalfSpace{N}, Hyperplane{N}, Line{N}}}
ρ(::AbstractVector{N}, ::Intersection{N, S1, S2}) where {N<:Real, S1<:LazySet{N}, S2<:AbstractPolytope{N}}
isempty(::Intersection)
σ(::AbstractVector{N}, ::Intersection{N}) where {N<:Real}
∈(::AbstractVector{N}, ::Intersection{N}) where {N<:Real}
isempty_known(::Intersection)
set_isempty!(::Intersection, ::Bool)
swap(::Intersection)
use_precise_ρ
_line_search
_projection
```

Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{N}) where {N<:Real})

#### Intersection cache

```@docs
IntersectionCache
```

### ``n``-ary Intersection

```@docs
IntersectionArray
dim(::IntersectionArray)
σ(::AbstractVector{N}, ::IntersectionArray{N}) where {N<:Real}
∈(::AbstractVector{N}, ::IntersectionArray{N}) where {N<:Real}
array(::IntersectionArray{N, S}) where {N<:Real, S<:LazySet{N}}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{N}) where {N<:Real})

## Minkowski Sum

### Binary Minkowski Sum

```@docs
MinkowskiSum
⊕(::LazySet, ::LazySet)
+(::LazySet, ::LazySet)
dim(::MinkowskiSum)
ρ(::AbstractVector{N}, ::MinkowskiSum{N}) where {N<:Real}
σ(::AbstractVector{N}, ::MinkowskiSum{N}) where {N<:Real}
isempty(::MinkowskiSum)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{N}) where {N<:Real})

### ``n``-ary Minkowski Sum

```@docs
MinkowskiSumArray
dim(::MinkowskiSumArray)
ρ(::AbstractVector{N}, ::MinkowskiSumArray{N}) where {N<:Real}
σ(::AbstractVector{N}, ::MinkowskiSumArray{N}) where {N<:Real}
isempty(::MinkowskiSumArray)
array(::MinkowskiSumArray{N, S}) where {N<:Real, S<:LazySet{N}}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{N}) where {N<:Real})

### ``n``-ary Minkowski Sum with cache

```@docs
CacheMinkowskiSum
dim(::CacheMinkowskiSum)
σ(::AbstractVector{N}, ::CacheMinkowskiSum{N}) where {N<:Real}
isempty(::CacheMinkowskiSum)
array(::CacheMinkowskiSum{N, S}) where {N<:Real, S<:LazySet{N}}
forget_sets!(::CacheMinkowskiSum)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{N}) where {N<:Real})

## Maps

### Linear Map

```@docs
LinearMap
*(::AbstractMatrix{N}, ::LazySet{N}) where {N<:Real}
*(::N, ::LazySet{N}) where {N<:Real}
*(::N, ::LM) where {N<:Real, LM<:LinearMap{N}}
*(::AbstractMatrix{N}, ::ZeroSet{N}) where {N<:Real}
dim(::LinearMap)
ρ(::AbstractVector{N}, ::LinearMap{N}) where {N<:Real}
σ(::AbstractVector{N}, ::LinearMap{N}) where {N<:Real}
∈(::AbstractVector{N}, ::LinearMap{N}) where {N<:Real}
an_element(::LinearMap{N}) where {N<:Real}
isempty(::LinearMap)
vertices_list(::LinearMap{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

### Exponential Map

```@docs
ExponentialMap
dim(::ExponentialMap)
ρ(::AbstractVector{N}, ::ExponentialMap{N}) where {N<:Real}
σ(::AbstractVector{N}, ::ExponentialMap{N}) where {N<:Real}
∈(::AbstractVector{N}, ::ExponentialMap{N}) where {N<:Real}
isempty(::ExponentialMap)
vertices_list(::ExponentialMap{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{N}) where {N<:Real})

```@docs
ExponentialProjectionMap
dim(::ExponentialProjectionMap)
σ(::AbstractVector{N}, ::ExponentialProjectionMap{N}) where {N<:Real}
isempty(::ExponentialProjectionMap)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{N}) where {N<:Real})

```@docs
SparseMatrixExp
*(::SparseMatrixExp{N}, ::LazySet{N}) where {N<:Real}
get_row(::SparseMatrixExp, ::Int)
```

```@docs
ProjectionSparseMatrixExp
*(::ProjectionSparseMatrixExp, ::LazySet)
```

## Symmetric Interval Hull

```@docs
SymmetricIntervalHull
dim(::SymmetricIntervalHull)
σ(::AbstractVector{N}, ::SymmetricIntervalHull{N}) where {N<:Real}
center(::SymmetricIntervalHull{N}) where {N<:Real}
radius_hyperrectangle(::SymmetricIntervalHull{N}) where {N<:Real}
radius_hyperrectangle(::SymmetricIntervalHull{N}, ::Int) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`singleton_list`](@ref singleton_list(::AbstractPolytope{N}) where {N<:Real})
* [`linear_map`](@ref linear_map(::AbstractMatrix{N}, ::AbstractPolytope{N}) where {N<:Real})

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope{N}) where {N<:Real})

Inherited from [`AbstractHyperrectangle`](@ref):
* [`∈`](@ref ∈(::AbstractVector{N}, ::AbstractHyperrectangle{N}) where {N<:Real})
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real))
* [`vertices_list`](@ref vertices_list(::AbstractHyperrectangle{N}) where {N<:Real})
* [`high`](@ref high(::AbstractHyperrectangle{N}) where {N<:Real})
* [`low`](@ref low(::AbstractHyperrectangle{N}) where {N<:Real})
