"""
    complement(∅::EmptySet{N}) where {N}

Return the complement of an empty set.

### Input

- `∅` -- empty set

### Output

The universe of the same dimension.
"""
function complement(∅::EmptySet{N}) where {N}
    return Universe{N}(dim(∅))
end
