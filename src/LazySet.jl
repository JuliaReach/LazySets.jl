import Base.LinAlg:norm

export LazySet,
       ρ, support_function,
       σ, support_vector,
       dim,
       norm,
       radius,
       diameter,
       an_element,
       neutral,
       absorbing

"""
    LazySet{N}

Abstract type for set representations and lazy operations between sets.
Most subtypes of `LazySet` correspond to convex sets; however, convexity is not
a restriction for subtyping `LazySet`.

### Notes

`LazySet` types should be parameterized with a type `N`, typically
`N<:Real`, for using different numeric types.

Every concrete `LazySet` must define the following functions:
- `σ(d::AbstractVector{N}, S::LazySet)::AbstractVector{N}` -- the
    support vector of `S` in a given direction `d`
- `dim(S::LazySet)::Int` -- the ambient dimension of `S`

```jldoctest
julia> subtypes(LazySet)
17-element Array{Union{DataType, UnionAll},1}:
 LazySets.AbstractPointSymmetric
 LazySets.AbstractPolytope
 LazySets.CartesianProduct
 LazySets.CartesianProductArray
 LazySets.ConvexHull
 LazySets.ConvexHullArray
 LazySets.EmptySet
 LazySets.ExponentialMap
 LazySets.ExponentialProjectionMap
 LazySets.HalfSpace
 LazySets.Hyperplane
 LazySets.Intersection
 LazySets.Line
 LazySets.LinearMap
 LazySets.MinkowskiSum
 LazySets.MinkowskiSumArray
 LazySets.PolynomialZonotope
```
"""
abstract type LazySet{N} end


# --- common LazySet functions ---


"""
    ρ(d::AbstractVector{N}, S::LazySet{N})::N where {N<:Real}

Evaluate the support function of a set in a given direction.

### Input

- `d` -- direction
- `S` -- set

### Output

The support function of the set `S` for the direction `d`.
"""
function ρ(d::AbstractVector{N}, S::LazySet{N})::N where {N<:Real}
    return dot(d, σ(d, S))
end

"""
    support_function

Alias for the support function ρ.
"""
const support_function = ρ

"""
    σ

Function to compute the support vector σ.
"""
function σ end

"""
    support_vector

Alias for the support vector σ.
"""
const support_vector = σ


"""
    isconvex(S::LazySet)::Bool

Sufficient check if a set is convex.

### Input

- `S` -- set

### Output

`false` by default.
"""
function isconvex(S::LazySet)::Bool
    return false
end


"""
    norm(S::LazySet, [p]::Real=Inf)

Return the norm of a set.
It is the norm of the enclosing ball (of the given ``p``-norm) of minimal volume
that is centered in the origin.

### Input

- `S` -- set
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the norm.
"""
function norm(S::LazySet, p::Real=Inf)
    if p == Inf
        return norm(Approximations.ballinf_approximation(S), p)
    else
        error("the norm for this value of p=$p is not implemented")
    end
end

"""
    radius(S::LazySet, [p]::Real=Inf)

Return the radius of a set.
It is the radius of the enclosing ball (of the given ``p``-norm) of minimal
volume with the same center.

### Input

- `S` -- set
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the radius.
"""
function radius(S::LazySet, p::Real=Inf)
    if p == Inf
        return radius(Approximations.ballinf_approximation(S)::BallInf, p)
    else
        error("the radius for this value of p=$p is not implemented")
    end
end

"""
    diameter(S::LazySet, [p]::Real=Inf)

Return the diameter of a set.
It is the maximum distance between any two elements of the set, or,
equivalently, the diameter of the enclosing ball (of the given ``p``-norm) of
minimal volume with the same center.

### Input

- `S` -- set
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the diameter.
"""
function diameter(S::LazySet, p::Real=Inf)
    if p == Inf
        return radius(S, p) * 2
    else
        error("the diameter for this value of p=$p is not implemented")
    end
end


"""
    an_element(S::LazySet{N})::AbstractVector{N} where {N<:Real}

Return some element of a set.

### Input

- `S` -- set

### Output

An element of a set.
"""
function an_element(S::LazySet{N})::AbstractVector{N} where {N<:Real}
    return σ(sparsevec([1], [one(N)], dim(S)), S)
end
