function σ(d::AbstractVector, L::Line2D)
    v, unbounded = _σ_hyperplane_halfspace(d, L.a, L.b; error_unbounded=true,
                                           halfspace=false)
    return v
end
