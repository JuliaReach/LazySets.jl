import Base: extrema

export LazySet,
       low,
       high

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
9-element Vector{Any}:
 AbstractPolynomialZonotope
 Complement
 ConvexSet
 LazySets.AbstractStar
 QuadraticMap
 Rectification
 UnionSet
 UnionSetArray
 VPolygonNC
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
VPolygonNC
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

# Note: this method cannot be documented due to a bug in Julia
function low(X::LazySet, i::Int)
    return _low(X, i)
end

function _low(X::LazySet{N}, i::Int) where {N}
    n = dim(X)
    d = SingleEntryVector(i, n, -one(N))
    return -ρ(d, X)
end

"""
    low(X::LazySet)

Return a vector with the lowest coordinates of the set in each canonical
direction.

### Input

- `X` -- set

### Output

A vector with the lower coordinate of the set in each dimension.

### Notes

See also [`low(X::LazySet, i::Int)`](@ref).

The result is the lowermost corner of the box approximation, so it is not
necessarily contained in `X`.
"""
function low(X::LazySet)
    n = dim(X)
    return [low(X, i) for i in 1:n]
end

# Note: this method cannot be documented due to a bug in Julia
function high(X::LazySet, i::Int)
    return _high(X, i)
end

function _high(X::LazySet{N}, i::Int) where {N}
    n = dim(X)
    d = SingleEntryVector(i, n, one(N))
    return ρ(d, X)
end

"""
    high(X::LazySet)

Return a vector with the highest coordinate of the set in each canonical
direction.

### Input

- `X` -- set

### Output

A vector with the highest coordinate of the set in each dimension.

### Notes

See also [`high(X::LazySet, i::Int)`](@ref).

The result is the uppermost corner of the box approximation, so it is not
necessarily contained in `X`.
"""
function high(X::LazySet)
    n = dim(X)
    return [high(X, i) for i in 1:n]
end

"""
    extrema(X::LazySet, i::Int)

Return the lower and higher coordinate of a set in a given dimension.

### Input

- `X` -- set
- `i` -- dimension of interest

### Output

The lower and higher coordinate of the set in the given dimension.

### Notes

The result is equivalent to `(low(X, i), high(X, i))`, but sometimes it can be
computed more efficiently.

### Algorithm

The bounds are computed with `low` and `high`.
"""
function extrema(X::LazySet, i::Int)
    l = low(X, i)
    h = high(X, i)
    return (l, h)
end

"""
    extrema(X::LazySet)

Return two vectors with the lowest and highest coordinate of `X` in each
canonical direction.

### Input

- `X` -- set

### Output

Two vectors with the lowest and highest coordinates of `X` in each dimension.

### Notes

The result is equivalent to `(low(X), high(X))`, but sometimes it can be
computed more efficiently.

The resulting points are the lowermost and uppermost corners of the box
approximation, so they are not necessarily contained in `X`.

### Algorithm

The bounds are computed with `low` and `high`.
"""
function extrema(X::LazySet)
    l = low(X)
    h = high(X)
    return (l, h)
end

"""
    plot_vlist(X::S, ε::Real) where {S<:LazySet}

Return a list of vertices used for plotting a two-dimensional set.

### Input

- `X` -- two-dimensional set
- `ε` -- precision parameter

### Output

A list of vertices of a polygon `P`.
For convex `X`, `P` usually satisfies that the Hausdorff distance to `X` is less
than `ε`.
"""
function plot_vlist(X::S, ε::Real) where {S<:LazySet}
    @assert isconvextype(S) "can only plot convex sets"

    P = overapproximate(X, ε)
    return convex_hull(vertices_list(P))
end
