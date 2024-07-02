"""
    radius_hyperrectangle(B::BallInf, i::Int)

Return the box radius of a ball in the infinity norm in a given dimension.

### Input

- `B` -- ball in the infinity norm
- `i` -- dimension of interest

### Output

The box radius of the ball in the infinity norm in the given dimension.
"""
function radius_hyperrectangle(B::BallInf, i::Int)
    @assert 1 <= i <= dim(B) "cannot compute the radius of a " *
                             "$(dim(B))-dimensional set in dimension $i"
    return B.radius
end

"""
    radius_hyperrectangle(B::BallInf)

Return the box radius of a ball in the infinity norm.

### Input

- `B` -- ball in the infinity norm

### Output

The box radius of the ball in the infinity norm, which is the same in every
dimension.
"""
function radius_hyperrectangle(B::BallInf)
    return fill(B.radius, dim(B))
end
