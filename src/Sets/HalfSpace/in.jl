"""
# Extended help

    ∈(x::AbstractVector, hs::HalfSpace)

### Algorithm

This implementation checks whether ``x`` satisfies ``a⋅x ≤ b``.
"""
function ∈(x::AbstractVector, hs::HalfSpace)
    return _leq(dot(x, hs.a), hs.b)
end
