"""
    distance(x::AbstractVector, H::Hyperplane)

Compute the distance between point `x` and hyperplane `H` with respect to the
Euclidean norm.

### Input

- `x` -- vector
- `H` -- hyperplane

### Output

A scalar representing the distance between point `x` and hyperplane `H`.
"""
@commutative function distance(x::AbstractVector, H::Hyperplane)
    N = promotetype(eltype(x), eltype(H))
    a, b = _normalize_halfspace(H, N(2))
    return abs(dot(x, a) - b)
end
