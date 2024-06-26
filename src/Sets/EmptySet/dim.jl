"""
    dim(∅::EmptySet)

Return the dimension of an empty set.

### Input

- `∅` -- an empty set

### Output

The dimension of the empty set.
"""
function dim(∅::EmptySet)
    return ∅.dim
end
