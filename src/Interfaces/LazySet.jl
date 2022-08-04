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
