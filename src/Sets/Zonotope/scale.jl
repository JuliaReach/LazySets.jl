"""
    scale!(α::Real, Z::Zonotope)

Concrete scaling of a zonotope modifying `Z` in-place.

### Input

- `α` -- scalar
- `Z` -- zonotope

### Output

The zonotope obtained by applying the numerical scale to the center and
generators of ``Z``.
"""
function scale!(α::Real, Z::Zonotope)
    Z.center .*= α
    Z.generators .*= α
    return Z
end
