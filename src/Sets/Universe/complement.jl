"""
    complement(U::Universe{N}) where {N}

Return the complement of an universe.

### Input

- `∅` -- universe

### Output

The empty set of the same dimension.
"""
function complement(U::Universe{N}) where {N}
    require(@__MODULE__, :LazySets; fun_name="complement")

    return EmptySet{N}(dim(U))
end
