"""
    low(∅::EmptySet)

Return a vector with the lowest coordinates of an empty set in each canonical
direction.

### Input

- `∅` -- empty set

### Output

An error.

### Notes

See also [`low(∅::EmptySet, i::Int)`](@ref).
"""
function low(::EmptySet)
    return error("the lower bound of an empty set is undefined")
end

"""
    low(∅::EmptySet, i::Int)

Return the lowest coordinate of an empty set in the given direction.

### Input

- `∅` -- empty set
- `i` -- dimension of interest

### Output

An error.
"""
function low(::EmptySet, ::Int)
    return error("the lower bound of an empty set is undefined")
end
