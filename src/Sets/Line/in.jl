"""
# Extended help

    ∈(x::AbstractVector, L::Line)

### Algorithm

The point ``x`` belongs to the line ``L : p + λd`` if and only if ``x - p`` is
proportional to the direction ``d``.
"""
@validate function ∈(x::AbstractVector, L::Line)
    # the following check is necessary because the zero vector is a special case
    _isapprox(x, L.p) && return true

    return first(ismultiple(x - L.p, L.d))
end
