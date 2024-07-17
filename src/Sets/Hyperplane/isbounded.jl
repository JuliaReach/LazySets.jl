"""
    isbounded(H::Hyperplane)

Check whether a hyperplane is bounded.

### Input

- `H` -- hyperplane

### Output

`true` iff `H` is one-dimensional.
"""
function isbounded(H::Hyperplane)
    return dim(H) == 1
end
