"""
# Extended help

    radius(B::BallInf, [p]::Real=Inf)

### Notes

The result is defined as the radius of the enclosing ball of the given
``p``-norm of minimal volume with the same center.
"""
@validate function radius(B::BallInf, p::Real=Inf)
    return (p == Inf) ? B.radius : norm(fill(B.radius, dim(B)), p)
end
