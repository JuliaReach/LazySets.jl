"""
    dim(cup::UnionSet)

Return the dimension of the union of two sets.

### Input

- `cup` -- union of two sets

### Output

The ambient dimension of the union of two sets.
"""
function dim(cup::UnionSet)
    return dim(cup.X)
end
