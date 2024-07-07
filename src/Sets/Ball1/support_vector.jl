"""
    σ(d::AbstractVector, B::Ball1)

Return the support vector of a ball in the 1-norm in the given direction.

### Input

- `d` -- direction
- `B` -- ball in the 1-norm

### Output

The support vector in the given direction.
"""
function σ(d::AbstractVector, B::Ball1)
    res = copy(B.center)
    imax = argmaxabs(d)
    res[imax] += sign(d[imax]) * B.radius
    return res
end
