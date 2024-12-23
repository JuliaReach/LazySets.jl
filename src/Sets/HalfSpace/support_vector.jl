"""
# Extended help

    σ(d::AbstractVector, hs::HalfSpace)

### Output

The support vector in the given direction, which is only defined in the
following two cases:
1. The direction has norm zero.
2. The direction is (a multiple of) the normal direction of the half-space.
In both cases the result is any point on the boundary (the defining hyperplane).
Otherwise this function throws an error.
"""
function σ(d::AbstractVector, hs::HalfSpace)
    v, unbounded = _σ_hyperplane_halfspace(d, hs.a, hs.b; error_unbounded=true,
                                           halfspace=true)
    return v
end
