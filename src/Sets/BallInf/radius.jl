"""
    radius(B::BallInf, [p]::Real=Inf)

Return the radius of a ball in the infinity norm.

### Input

- `B` -- ball in the infinity norm
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the radius.

### Notes

The result is defined as the radius of the enclosing ball of the given
``p``-norm of minimal volume with the same center.
"""
function radius(B::BallInf, p::Real=Inf)
    return (p == Inf) ? B.radius : norm(fill(B.radius, dim(B)), p)
end
