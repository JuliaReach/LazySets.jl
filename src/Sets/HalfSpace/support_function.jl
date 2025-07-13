"""
# Extended help

    ρ(d::AbstractVector, hs::HalfSpace)

### Output

Unless the direction is (a multiple of) the normal direction of the half-space,
the result is `Inf`.
"""
@validate function ρ(d::AbstractVector, hs::HalfSpace)
    v, unbounded = _σ_hyperplane_halfspace(d, hs.a, hs.b; error_unbounded=false,
                                           halfspace=true)
    if unbounded
        N = promote_type(eltype(d), eltype(hs))
        return N(Inf)
    end
    return dot(d, v)
end
