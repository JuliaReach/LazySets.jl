"""
    constraints_list(U::Universe{N}) where {N}

Return the list of constraints defining a universe.

### Input

- `U` -- universe

### Output

The empty list of constraints, as the universe is unconstrained.
"""
function constraints_list(U::Universe{N}) where {N}
    return HalfSpace{N,Vector{N}}[]
end
