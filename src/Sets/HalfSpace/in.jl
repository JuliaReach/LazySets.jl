"""
# Extended help

    in(x::AbstractVector, hs::HalfSpace)

### Algorithm

This implementation checks whether ``x`` satisfies ``a⋅x ≤ b``.
"""
@validate function in(x::AbstractVector, hs::HalfSpace)
    return _leq(dot(x, hs.a), hs.b)
end
