"""
    constraints_list(lm::LinearMap)

Return the list of constraints of a (polyhedral) linear map.

### Input

- `lm` -- linear map

### Output

The list of constraints of the linear map.

### Notes

We assume that the underlying set `X` is polyhedral, i.e., offers a method
`constraints_list(X)`.

### Algorithm

We fall back to a concrete set representation by applying `linear_map`.
"""
function constraints_list(lm::LinearMap)
    return constraints_list(linear_map(lm.M, lm.X))
end
