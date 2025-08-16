"""
    radius_hyperrectangle(H::Hyperrectangle)

Return the box radius of a hyperrectangle in every dimension.

### Input

- `H` -- hyperrectangle

### Output

The box radius of the hyperrectangle.
"""
function radius_hyperrectangle(H::Hyperrectangle)
    return H.radius
end

"""
    radius_hyperrectangle(H::Hyperrectangle, i::Int)

Return the box radius of a hyperrectangle in a given dimension.

### Input

- `H` -- hyperrectangle
- `i` -- dimension of interest

### Output

The box radius in the given dimension.
"""
@validate function radius_hyperrectangle(H::Hyperrectangle, i::Int)
    return H.radius[i]
end
