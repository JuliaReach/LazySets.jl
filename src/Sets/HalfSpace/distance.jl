"""
    distance(x::AbstractVector, H::HalfSpace)

Compute the distance between point `x` and half-space `H` with respect to the
Euclidean norm.

### Input

- `x` -- vector
- `H` -- half-space

### Output

A scalar representing the distance between point `x` and half-space `H`.
"""
@commutative function distance(x::AbstractVector, H::HalfSpace)
    N = promote_type(eltype(x), eltype(H))
    a, b = _normalize_halfspace(H, N(2))
    return max(dot(x, a) - b, zero(N))
end
