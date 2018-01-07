# =============================
# Support function of a LazySet
# =============================

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
support_function = ρ

"""
    support_vector

Alias for the support vector σ.
"""
support_vector = σ
