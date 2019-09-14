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

## Cartesian Product

### Binary Cartesian Product

```@docs
CartesianProduct
×(::LazySet, ::LazySet)
*(::LazySet, ::LazySet)
swap(::CartesianProduct)
dim(::CartesianProduct)
ρ(::AbstractVector{N}, ::CartesianProduct{N}) where {N<:Real}
σ(::AbstractVector{N}, ::CartesianProduct{N}) where {N<:Real}
isbounded(::CartesianProduct)
∈(::AbstractVector{N}, ::CartesianProduct{N}) where {N<:Real}
isempty(::CartesianProduct)
constraints_list(::CartesianProduct{N}) where {N<:Real}
vertices_list(::CartesianProduct{N}) where {N<:Real}
linear_map(M::AbstractMatrix{N}, cp::CartesianProduct{N}) where {N<:Real}
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
isbounded(::CartesianProductArray)
∈(::AbstractVector{N}, ::CartesianProductArray{N}) where {N<:Real}
isempty(::CartesianProductArray)
constraints_list(::CartesianProductArray{N}) where {N<:Real}
vertices_list(::CartesianProductArray{N}) where {N<:Real}
linear_map(M::AbstractMatrix{N}, cpa::CartesianProductArray{N}) where {N<:Real}
array(::CartesianProductArray{N, S}) where {N<:Real, S<:LazySet{N}}
block_structure(cpa::CartesianProductArray{N}) where {N}
block_to_dimension_indices(cpa::CartesianProductArray{N}, vars::Vector{Int}) where {N}
substitute_blocks(low_dim_cpa::CartesianProductArray{N}, orig_cpa::CartesianProductArray{N},
blocks::Vector{Tuple{Int,Int}}) where {N}
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
swap(::ConvexHull)
dim(::ConvexHull)
ρ(::AbstractVector{N}, ::ConvexHull{N}) where {N<:Real}
σ(::AbstractVector{N}, ::ConvexHull{N}) where {N<:Real}
isbounded(::ConvexHull)
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
isbounded(::ConvexHullArray)
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
convex_hull(::Vector{VN}) where {N<:Real, VN<:AbstractVector{N}}
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
ρ(::AbstractVector{N}, ::Intersection{N, S1, S2}) where {N<:Real, S1<:LazySet{N}, S2<:AbstractPolyhedron{N}}
ρ(::AbstractVector{N}, ::Intersection{N, S1, S2}) where {N<:Real, S1<:AbstractPolyhedron{N}, S2<:AbstractPolyhedron{N}}
σ(::AbstractVector{N}, ::Intersection{N}) where {N<:Real}
isbounded(::Intersection)
isempty(::Intersection)
∈(::AbstractVector{N}, ::Intersection{N}) where {N<:Real}
constraints_list(::Intersection{N}) where {N<:Real}
isempty_known(::Intersection)
set_isempty!(::Intersection, ::Bool)
swap(::Intersection{N, S1, S2}) where {N<:Real, S1, S2}
use_precise_ρ
_line_search
_projection
linear_map(::AbstractMatrix{N}, ::Intersection{N}) where {N}
plot_recipe(::Intersection{N}, ::N=zero(N), ::Int=40) where {N<:Real}
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::Intersection{N}, ::N=zero(N), ::Int=40) where {N<:Real}
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
isbounded(::IntersectionArray)
∈(::AbstractVector{N}, ::IntersectionArray{N}) where {N<:Real}
array(::IntersectionArray{N, S}) where {N<:Real, S<:LazySet{N}}
constraints_list(::IntersectionArray{N}) where {N<:Real}
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
swap(::MinkowskiSum)
dim(::MinkowskiSum)
ρ(::AbstractVector{N}, ::MinkowskiSum{N}) where {N<:Real}
σ(::AbstractVector{N}, ::MinkowskiSum{N}) where {N<:Real}
isbounded(::MinkowskiSum)
isempty(::MinkowskiSum)
constraints_list(::MinkowskiSum)
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
isbounded(::MinkowskiSumArray)
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
isbounded(::CacheMinkowskiSum)
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
*(::Union{AbstractMatrix, UniformScaling, AbstractVector, Real}, ::LazySet)
dim(::LinearMap)
ρ(::AbstractVector{N}, ::LinearMap{N}) where {N<:Real}
σ(::AbstractVector{N}, ::LinearMap{N}) where {N<:Real}
∈(::AbstractVector{N}, ::LinearMap{N}) where {N<:Real}
an_element(::LinearMap{N}) where {N<:Real}
isbounded(::LinearMap)
isempty(::LinearMap)
vertices_list(::LinearMap{N}) where {N<:Real}
constraints_list(::LinearMap{N}) where {N<:Real}
linear_map(::AbstractMatrix{N}, ::LinearMap{N}) where {N<:Real}
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
isbounded(::ExponentialMap)
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
isbounded(::ExponentialProjectionMap)
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

### Reset Map

```@docs
ResetMap
dim(::ResetMap)
ρ(::AbstractVector{N}, ::ResetMap{N}) where {N<:Real}
σ(::AbstractVector{N}, ::ResetMap{N}) where {N<:Real}
an_element(::ResetMap)
isempty(::ResetMap)
get_A(::ResetMap{N}) where {N<:Real}
get_b(::ResetMap{N}) where {N<:Real}
constraints_list(::ResetMap{N}) where {N<:Real}
constraints_list(::ResetMap{N, S}) where {N<:Real, S<:AbstractHyperrectangle}
```

Inherited from [`LazySet`](@ref):
* [`isbounded`](@ref isbounded(::LinearMap))
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

### Translation

```@docs
Translation
+(X::LazySet, v::AbstractVector)
⊕(X::LazySet, v::AbstractVector)
dim(::Translation)
ρ(::AbstractVector{N}, ::Translation{N}) where {N<:Real}
σ(::AbstractVector{N}, ::Translation{N}) where {N<:Real}
an_element(::Translation)
isempty(::Translation)
constraints_list(::Translation{N}) where {N<:Real}
linear_map(::AbstractMatrix{N}, ::Translation{N}) where {N<:Real}
∈(::AbstractVector{N}, ::Translation{N}) where {N<:Real}
```

### Affine Map

```@docs
AffineMap
dim(::AffineMap)
σ(::AbstractVector{N}, ::AffineMap{N}) where {N<:Real}
ρ(::AbstractVector{N}, ::AffineMap{N}) where {N<:Real}
an_element(::AffineMap)
isempty(::AffineMap)
isbounded(::AffineMap)
∈(::AbstractVector{N}, ::AffineMap{N}) where {N<:Real}
vertices_list(::AffineMap{N}) where {N<:Real}
constraints_list(::AffineMap{N}) where {N<:Real}
linear_map(::AbstractMatrix{N}, ::AffineMap{N}) where {N<:Real}
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
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`singleton_list`](@ref singleton_list(::AbstractPolytope{N}) where {N<:Real})

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope{N}) where {N<:Real})

Inherited from [`AbstractZonotope`](@ref):
* [`ngens`](@ref ngens(::AbstractZonotope))
* [`linear_map`](@ref linear_map(::AbstractMatrix{N}, ::AbstractZonotope{N}) where {N<:Real})
* [`order`](@ref order(::AbstractZonotope))
* [`translate`](@ref translate(::AbstractZonotope{N}, ::AbstractVector{N}) where {N<:Real})

Inherited from [`AbstractHyperrectangle`](@ref):
* [`ρ`](@ref ρ(::AbstractVector{N}, ::AbstractHyperrectangle{N}) where {N<:Real})
* [`σ`](@ref σ(::AbstractVector{N}, ::AbstractHyperrectangle{N}) where {N<:Real})
* [`∈`](@ref ∈(::AbstractVector{N}, ::AbstractHyperrectangle{N}) where {N<:Real})
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real))
* [`constraints_list`](@ref constraints_list(::AbstractHyperrectangle{N}) where {N<:Real})
* [`vertices_list`](@ref vertices_list(::AbstractHyperrectangle{N}) where {N<:Real})
* [`high`](@ref high(::AbstractHyperrectangle{N}) where {N<:Real})
* [`low`](@ref low(::AbstractHyperrectangle{N}) where {N<:Real})
* [`generators`](@ref generators(::AbstractHyperrectangle))
* [`genmat`](@ref genmat(::AbstractHyperrectangle))

## Union

Note that the union of convex sets is generally not convex.
Hence these set types are not part of the convex-set family `LazySet`.

### Binary Set Union

```@docs
UnionSet
∪(::LazySet, ::LazySet)
swap(::UnionSet)
dim(::UnionSet)
σ(::AbstractVector{N}, ::UnionSet{N}; algorithm="support_vector") where {N<:Real}
ρ(::AbstractVector{N}, ::UnionSet{N}) where {N<:Real}
an_element(::UnionSet{N}) where {N<:Real}
∈(::AbstractVector{N}, ::UnionSet{N}) where {N<:Real}
isempty(::UnionSet)
isbounded(::UnionSet)
```

### ``n``-ary Set Union

```@docs
UnionSetArray
dim(::UnionSetArray)
array(::UnionSetArray{N, S}) where {N<:Real, S<:LazySet{N}}
σ(::AbstractVector{N}, ::UnionSetArray{N}; algorithm="support_vector") where {N<:Real}
ρ(::AbstractVector{N}, ::UnionSetArray{N}) where {N<:Real}
an_element(::UnionSetArray{N}) where {N<:Real}
∈(::AbstractVector{N}, ::UnionSetArray{N}) where {N<:Real}
isempty(::UnionSetArray)
isbounded(::UnionSetArray)
```

## Complement

Note that the complement of a convex set is generally not convex.
Hence this set type is not part of the convex-set family `LazySet`.

```@docs
Complement
dim(::Complement)
∈(::AbstractVector{N}, ::Complement{N}) where {N<:Real}
isempty(::Complement)
```

## Rectification

Note that the rectification of a convex set is generally not convex.
Hence this set type is not part of the convex-set family `LazySet`.

```@docs
Rectification
dim(::Rectification)
σ(::AbstractVector{N}, ::Rectification{N}) where {N<:Real}
σ(::AbstractVector{N}, ::Rectification{N, <:AbstractHyperrectangle{N}}) where {N<:Real}
σ(::AbstractVector{N}, ::Rectification{N, <:CartesianProduct{N}}) where {N<:Real}
σ(::AbstractVector{N}, ::Rectification{N, <:CartesianProductArray{N}}) where {N<:Real}
ρ(::AbstractVector{N}, ::Rectification{N}) where {N<:Real}
an_element(::Rectification{N}) where {N<:Real}
∈(::AbstractVector{N}, ::Rectification{N}) where {N<:Real}
isempty(::Rectification)
isbounded(::Rectification{N}) where {N<:Real}
to_union_of_projections(::Rectification{N}, ::Bool=false) where {N<:Real}
```

#### Rectification cache

```@docs
RectificationCache
```

## Interior

```@docs
is_interior_point(::AbstractVector, ::LazySet; p=Inf, ε=_TOL_F64.rtol)
```
