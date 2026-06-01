"""
    constraints_list(cp::CartesianProduct)

Return the list of constraints of a (polyhedral) Cartesian product.

### Input

- `cp` -- polyhedral Cartesian product

### Output

A list of constraints.
"""
function constraints_list(cp::CartesianProduct)
    return _constraints_list_cartesian_product(cp)
end
