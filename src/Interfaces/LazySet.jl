export LazySet

"""
    LazySet{N}

Abstract type for the set types in LazySets.

### Notes

`LazySet` types should be parameterized with a type `N`, typically `N<:Real`,
for using different numeric types.

Every concrete `LazySet` must define the following functions:
- `dim(S::LazySet)` -- the ambient dimension of `S`

The subtypes of `LazySet` (including abstract interfaces):

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(LazySet, false)
8-element Vector{Any}:
 AbstractPolynomialZonotope
 Complement
 ConvexSet
 LazySets.AbstractStar
 QuadraticMap
 Rectification
 UnionSet
 UnionSetArray
```

If we only consider *concrete* subtypes, then:

```jldoctest; setup = :(using LazySets: subtypes)
julia> concrete_subtypes = subtypes(LazySet, true);

julia> length(concrete_subtypes)
53

julia> println.(concrete_subtypes);
AffineMap
Ball1
Ball2
BallInf
Ballp
Bloating
CachedMinkowskiSumArray
CartesianProduct
CartesianProductArray
Complement
ConvexHull
ConvexHullArray
DensePolynomialZonotope
Ellipsoid
EmptySet
ExponentialMap
ExponentialProjectionMap
HParallelotope
HPolygon
HPolygonOpt
HPolyhedron
HPolytope
HalfSpace
Hyperplane
Hyperrectangle
Intersection
IntersectionArray
Interval
InverseLinearMap
LazySets.AbstractStar
Line
Line2D
LineSegment
LinearMap
MinkowskiSum
MinkowskiSumArray
QuadraticMap
Rectification
ResetMap
RotatedHyperrectangle
SimpleSparsePolynomialZonotope
Singleton
SparsePolynomialZonotope
Star
SymmetricIntervalHull
Translation
UnionSet
UnionSetArray
Universe
VPolygon
VPolytope
ZeroSet
Zonotope

```
"""
abstract type LazySet{N} end

"""
    isconvextype(X::Type{<:LazySet})

Check whether the given `LazySet` type is convex.

### Input

- `X` -- subtype of `LazySet`

### Output

`true` if the given set type is guaranteed to be convex by using only type
information, and `false` otherwise.

### Notes

Since this operation only acts on types (not on values), it can return false
negatives, i.e. there may be instances where the set is convex, even though the
answer of this function is `false`. The examples below illustrate this point.

### Examples

A ball in the infinity norm is always convex, hence we get:

```jldoctest convex_types
julia> isconvextype(BallInf)
true
```

For instance, the union (`UnionSet`) of two sets may in general be either convex
or not, since convexity cannot be decided by just using type information.
Hence, `isconvextype` returns `false` if `X` is `Type{<:UnionSet}`.

```jldoctest convex_types
julia> isconvextype(UnionSet)
false
```

However, the type parameters from the set operations allow to decide convexity
in some cases, by falling back to the convexity of the type of its arguments.
Consider for instance the lazy intersection. The intersection of two convex sets
is always convex, hence we can get:

```jldoctest convex_types
julia> isconvextype(Intersection{Float64, BallInf{Float64}, BallInf{Float64}})
true
```
"""
isconvextype(X::Type{<:LazySet}) = false
