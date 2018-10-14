import Base.==

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
elements ``x, y ∈ S`` and ``0 ≤ λ ≤ 1`` it holds that ``λ·x + (1-λ)·y ∈ S``.

### Notes

`LazySet` types should be parameterized with a type `N`, typically `N<:Real`,
for using different numeric types.

Every concrete `LazySet` must define the following functions:
- `σ(d::AbstractVector{N}, S::LazySet{N}) where {N<:Real}` -- the support vector
    of `S` in a given direction `d`; note that the numeric type `N` of `d` and
    `S` must be identical; for some set types `N` may be more restrictive than
    `Real`
- `dim(S::LazySet)::Int` -- the ambient dimension of `S`

```jldoctest
julia> subtypes(LazySet)
19-element Array{Any,1}:
 AbstractCentrallySymmetric
 AbstractPolytope
 CacheMinkowskiSum
 CartesianProduct
 CartesianProductArray
 ConvexHull
 ConvexHullArray
 EmptySet
 ExponentialMap
 ExponentialProjectionMap
 HPolyhedron
 HalfSpace
 Hyperplane
 Intersection
 IntersectionArray
 Line
 LinearMap
 MinkowskiSum
 MinkowskiSumArray
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

### Notes

The numeric type of the direction and the set must be identical.
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
    return radius(S, p) * 2
end


"""
    an_element(S::LazySet{N}) where {N<:Real}

Return some element of a convex set.

### Input

- `S` -- convex set

### Output

An element of a convex set.
"""
function an_element(S::LazySet{N}) where {N<:Real}
    return σ(sparsevec([1], [one(N)], dim(S)), S)
end


"""
    ==(X::LazySet, Y::LazySet)

Return whether two LazySets of the same type are exactly equal by recursively
comparing their fields until a mismatch is found.

### Input

- `X` -- any `LazySet`
- `Y` -- another `LazySet` of the same type as `X`

### Output

- `true` iff `X` is equal to `Y`.

### Notes

The check is purely syntactic and the sets need to have the same base type.
I.e. `X::VPolytope == Y::HPolytope` returns `false` even if `X` and `Y` represent the
same polytope. However `X::HPolytope{Int64} == Y::HPolytope{Float64}` is a valid comparison.

### Examples
```jldoctest
julia> HalfSpace([1], 1) == HalfSpace([1], 1)
true

julia> HalfSpace([1], 1) == HalfSpace([1.0], 1.0)
true

julia> Ball1([0.], 1.) == Ball2([0.], 1.)
false
```
"""
function ==(X::LazySet, Y::LazySet)
    # if the common supertype of X and Y is abstract, they cannot be compared
    if Compat.isabstracttype(promote_type(typeof(X), typeof(Y)))
        return false
    end

    for f in fieldnames(typeof(X))
        if getfield(X, f) != getfield(Y, f)
            return false
        end
    end

    return true
end
