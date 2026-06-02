"""
    isempty(cp::CartesianProduct)

Check whether a Cartesian product is empty.

### Input

- `cp` -- Cartesian product

### Output

`true` iff any of the sub-blocks is empty.
"""
function isempty(cp::CartesianProduct)
    return isempty(cp.X) || isempty(cp.Y)
end
