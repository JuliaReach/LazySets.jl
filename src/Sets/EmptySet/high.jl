"""
    high(∅::EmptySet)

Return a vector with the highest coordinates of an empty set in each canonical
direction.

### Input

- `∅` -- empty set

### Output

An error.

### Notes

See also [`high(∅::EmptySet, i::Int)`](@ref).
"""
function high(::EmptySet)
    return error("the upper bound of an empty set is undefined")
end

"""
    high(∅::EmptySet, i::Int)

Return the highest coordinate of an empty set in the given direction.

### Input

- `∅` -- empty set
- `i` -- dimension of interest

### Output

An error.
"""
function high(::EmptySet, ::Int)
    return error("the upper bound of an empty set is undefined")
end
