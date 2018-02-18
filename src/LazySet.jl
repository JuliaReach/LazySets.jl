import Base.LinAlg:norm,
       Base.∈

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

Abstract type for convex sets, i.e., sets characterized by a (possibly infinite)
intersection of halfspaces, or equivalently, sets ``S`` such that for any two
elements ``x, y ∈ S`` and ``0 ≤ λ ≤ 1`` it holds that ``λ x + (1-λ) y ∈ S``.

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
- `S` -- convex set

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
    norm(S::LazySet, [p]::Real=Inf)

Return the norm of a convex set.
It is the norm of the enclosing ball (of the given ``p``-norm) of minimal volume
that is centered in the origin.

### Input

- `S` -- convex set
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

Return the radius of a convex set.
It is the radius of the enclosing ball (of the given ``p``-norm) of minimal
volume with the same center.

### Input

- `S` -- convex set
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

Return the diameter of a convex set.
It is the maximum distance between any two elements of the set, or,
equivalently, the diameter of the enclosing ball (of the given ``p``-norm) of
minimal volume with the same center.

### Input

- `S` -- convex set
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

Return some element of a convex set.

### Input

- `S` -- convex set

### Output

An element of a convex set.
"""
function an_element(S::LazySet{N})::AbstractVector{N} where {N<:Real}
    return σ(sparsevec([1], [one(N)], dim(S)), S)
end

"""
    ∈(x::AbstractVector{N}, S::LazySet{N}, tolerance::N)::Bool
    where {N<:Real}

Check whether a given point is contained in a set with a given tolerance.

### Input

- `x` -- point/vector
- `S` -- set
- `tolerance` -- tolerance for when a point is still considered inside the set

### Output

`true` iff ``x ∈ S'``, where ``S'`` is the set ``S`` bloated by `tolerance`.

### Notes

This is a default implementation that tries different strategies which rely on
other implementations.

### Algorithm

First we try to check if ``X ∩ S ≠ ∅``, where ``X`` is a ball in the Euclidean
norm around ``x`` with radius `tolerance`.
If this fails (because the function is not implemented), we check containment in
the bloated set ``S' = S + B``, which is a Minkowski sum of ``S`` with a ball in
the Euclidean norm ``B`` around the origin with radius `tolerance`.
If this check also fails, we give up.
"""
function ∈(x::AbstractVector{N},
           S::LazySet{N},
           tolerance::N)::Bool where {N<:Real}
    @assert length(x) == dim(S)
    @assert tolerance >= 0
    try
        return !is_intersection_empty(Ball2(x, tolerance), S)
    catch MethodError
    end
    try
        return ∈(x, Ball2(zeros(N, length(x)), tolerance) + S)
    catch MethodError
    end
    error("containment check in $(typeof(S)) with tolerance is not implemented")
end
