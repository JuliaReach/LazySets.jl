"""
    normalize(hs::HalfSpace{N}, p::Real=N(2)) where {N}

Normalize a half-space.

### Input

- `hs` -- half-space
- `p`  -- (optional, default: `2`) norm

### Output

A new half-space whose normal direction ``a`` is normalized, i.e., such that
``‖a‖_p = 1`` holds.
"""
function normalize(hs::HalfSpace{N}, p::Real=N(2)) where {N}
    a, b = _normalize_halfspace(hs, p)
    return HalfSpace(a, b)
end

function _normalize_halfspace(H, p=2)
    nₐ = norm(H.a, p)
    a = LinearAlgebra.normalize(H.a, p)
    b = H.b / nₐ
    return a, b
end
