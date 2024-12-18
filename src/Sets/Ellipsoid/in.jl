"""
# Extended help

    ∈(x::AbstractVector, E::Ellipsoid)

### Algorithm

The point ``x`` belongs to the ellipsoid of center ``c`` and shape matrix ``Q``
if and only if

```math
(x-c)^\\mathrm{T} Q^{-1} (x-c) ≤ 1.
```
"""
function ∈(x::AbstractVector, E::Ellipsoid)
    @assert length(x) == dim(E) "cannot check membership of a vector of " *
                                "length $(length(x)) in an ellipsoid of dimension $(dim(E))"
    w = x - E.center
    Q = E.shape_matrix
    return dot(w, Q \ w) ≤ 1
end
