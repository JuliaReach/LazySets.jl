"""
# Extended help

    ∈(x::AbstractVector, L::Line)

### Algorithm

The point ``x`` belongs to the line ``L : p + λd`` if and only if ``x - p`` is
proportional to the direction ``d``.
"""
function ∈(x::AbstractVector, L::Line)
    @assert length(x) == dim(L) "expected the point and the line to have the " *
                                "same dimension, but they are $(length(x)) and $(dim(L)) respectively"

    # the following check is necessary because the zero vector is a special case
    _isapprox(x, L.p) && return true

    return first(ismultiple(x - L.p, L.d))
end
