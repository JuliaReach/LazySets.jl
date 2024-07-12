"""
    σ(d::AbstractVector, L::Line2D)

Return the support vector of a 2D line in a given direction.

### Input

- `d` -- direction
- `L` -- 2D line

### Output

The support vector in the given direction, which is defined the same way as for
the more general `Hyperplane`.
"""
function σ(d::AbstractVector, L::Line2D)
    v, unbounded = _σ_hyperplane_halfspace(d, L.a, L.b; error_unbounded=true,
                                           halfspace=false)
    return v
end
