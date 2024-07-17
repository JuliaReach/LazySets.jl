"""
    dim(H::Hyperplane)

Return the dimension of a hyperplane.

### Input

- `H` -- hyperplane

### Output

The ambient dimension of the hyperplane.
"""
function dim(H::Hyperplane)
    return length(H.a)
end
