# =============================
# Support function of a LazySet
# =============================

"""
    ρ(d::Vector{N}, S::LazySet)::N where {N<:Real}

Evaluate the support function of a set in a given direction.

### Input

- `d` -- a real vector, the direction investigated
- `S` -- a convex set

### Output

The support function of the set `S` for the direction `d`.
"""
function ρ(d::Vector{N}, S::LazySet)::N where {N<:Real}
    return dot(d, σ(d, S))
end

# aliases
support_function = ρ
support_vector = σ
