"""
    scale!(α::Real, P::DensePolynomialZonotope)

Scale a polynomial zonotope modified by a scale factor in-place.

### Input

- `α` -- scaling factor
- `P` -- polynomial zonotope

## Output

The modified polynomial zonotope.

### Algorithm

The center and generators are scaled by ``α``.
"""
function scale!(α::Real, P::DensePolynomialZonotope)
    P.c .*= α
    P.E .*= α
    P.F .*= α
    P.G .*= α
    return P
end
