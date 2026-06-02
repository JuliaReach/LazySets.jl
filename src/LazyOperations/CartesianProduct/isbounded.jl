"""
    isbounded(cp::CartesianProduct)

Check whether a Cartesian product is bounded.

### Input

- `cp` -- Cartesian product

### Output

`true` iff both wrapped sets are bounded.
"""
function isbounded(cp::CartesianProduct)
    return isbounded(cp.X) && isbounded(cp.Y)
end
