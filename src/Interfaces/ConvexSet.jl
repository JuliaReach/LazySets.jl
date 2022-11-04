export ConvexSet

"""
    ConvexSet{N} <: LazySet{N}

Abstract type for convex sets, i.e., sets characterized by a (possibly infinite)
intersection of halfspaces, or equivalently, sets ``S`` such that for any two
elements ``x, y ∈ S`` and ``0 ≤ λ ≤ 1`` it holds that ``λ·x + (1-λ)·y ∈ S``.
"""
abstract type ConvexSet{N} <: LazySet{N} end

isconvextype(X::Type{<:ConvexSet}) = true

"""
    low(X::ConvexSet{N}, i::Int) where {N}

Return the lower coordinate of a convex set in a given dimension.

### Input

- `X` -- convex set
- `i` -- dimension of interest

### Output

The lower coordinate of the set in the given dimension.
"""
function low(X::ConvexSet{N}, i::Int) where {N}
    # Note: this method is needed for documentation reasons
    # (see the method for LazySet)
    return _low(X, i)
end

"""
    high(X::ConvexSet{N}, i::Int) where {N}

Return the higher coordinate of a convex set in a given dimension.

### Input

- `X` -- convex set
- `i` -- dimension of interest

### Output

The higher coordinate of the set in the given dimension.
"""
function high(X::ConvexSet{N}, i::Int) where {N}
    # Note: this method is needed for documentation reasons
    # (see the method for LazySet)
    return _high(X, i)
end

# Note: this method is needed for documentation reasons
# (see the method for LazySet)
function an_element(S::ConvexSet{N}) where {N}
    return _an_element_lazySet(S)
end
