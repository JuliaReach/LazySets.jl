"""
# Extended help

    ∈(x::AbstractVector, L::Line2D)

### Algorithm

The point ``x`` belongs to the line if and only if ``a⋅x = b`` holds.
"""
@validate function ∈(x::AbstractVector, L::Line2D)
    return _isapprox(dot(L.a, x), L.b)
end
