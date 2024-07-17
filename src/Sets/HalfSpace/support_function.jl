"""
    ρ(d::AbstractVector, hs::HalfSpace)

Evaluate the support function of a half-space in a given direction.

### Input

- `d`  -- direction
- `hs` -- half-space

### Output

The support function of the half-space.
Unless the direction is (a multiple of) the normal direction of the half-space,
the result is `Inf`.
"""
function ρ(d::AbstractVector, hs::HalfSpace)
    v, unbounded = _σ_hyperplane_halfspace(d, hs.a, hs.b; error_unbounded=false,
                                           halfspace=true)
    if unbounded
        N = promote_type(eltype(d), eltype(hs))
        return N(Inf)
    end
    return dot(d, v)
end
