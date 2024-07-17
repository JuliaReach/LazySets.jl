"""
    normalize(H::Hyperplane{N}, p::Real=N(2)) where {N}

Normalize a hyperplane.

### Input

- `H` -- hyperplane
- `p` -- (optional, default: `2`) norm

### Output

A new hyperplane whose normal direction ``a`` is normalized, i.e., such that
``‖a‖_p = 1`` holds.
"""
function normalize(H::Hyperplane{N}, p::Real=N(2)) where {N}
    a, b = _normalize_halfspace(H, p)
    return Hyperplane(a, b)
end
