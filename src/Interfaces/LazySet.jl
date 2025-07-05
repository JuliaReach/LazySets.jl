export LazySet

"""
    LazySet{N}

Abstract type for the set types in LazySets.

### Notes

`LazySet` types should be parameterized with a type `N`, typically `N<:Real`,
for using different numeric types.

Every `LazySet` must implement the following function:

- `dim(X::LazySet)` -- the ambient dimension of `X`

While not strictly required, it is useful to implement the following function:

- `σ(d::AbstractVector, X::LazySet)` -- the support vector of `X` in a given
    direction `d`

Implementing the function

- `ρ(d::AbstractVector, X::LazySet)` -- the support function of `X` in a given
    direction `d`

is optional because there is a fallback implementation relying on `σ`.
However, for potentially unbounded sets (which includes most lazy set types)
this fallback cannot be used and an explicit implementation should be provided.

The subtypes of `LazySet` (including abstract interfaces):

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(LazySet, false)
17-element Vector{Any}:
 AbstractAffineMap
 AbstractPolynomialZonotope
 Bloating
 CachedMinkowskiSumArray
 CartesianProduct
 CartesianProductArray
 Complement
 ConvexSet
 Intersection
 IntersectionArray
 MinkowskiSum
 MinkowskiSumArray
 Polygon
 QuadraticMap
 Rectification
 UnionSet
 UnionSetArray
```

If we only consider *concrete* subtypes, then:

```jldoctest; setup = :(using LazySets: subtypes)
julia> concrete_subtypes = subtypes(LazySet, true);

julia> length(concrete_subtypes)
54

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
Line
Line2D
LineSegment
LinearMap
MinkowskiSum
MinkowskiSumArray
Polygon
QuadraticMap
Rectification
ResetMap
SimpleSparsePolynomialZonotope
Singleton
SparsePolynomialZonotope
Star
SymmetricIntervalHull
Tetrahedron
Translation
UnionSet
UnionSetArray
Universe
VPolygon
VPolytope
ZeroSet
Zonotope
ZonotopeMD
```
"""
abstract type LazySet{N} end
