"""
    ρ(d::AbstractVector, H::Hyperplane)

Evaluate the support function of a hyperplane in a given direction.

### Input

- `d` -- direction
- `H` -- hyperplane

### Output

The support function of the hyperplane.
If the set is unbounded in the given direction, the result is `Inf`.
"""
function ρ(d::AbstractVector, H::Hyperplane)
    v, unbounded = _σ_hyperplane_halfspace(d, H.a, H.b; error_unbounded=false,
                                           halfspace=false)
    if unbounded
        N = promote_type(eltype(d), eltype(H))
        return N(Inf)
    end
    return dot(d, v)
end
