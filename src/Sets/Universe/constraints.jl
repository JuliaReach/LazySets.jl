"""
    constraints(U::Universe{N}) where {N}

Construct an iterator over the constraints of a universe.

### Input

- `U` -- universe

### Output

The empty iterator, as the universe is unconstrained.
"""
function constraints(U::Universe{N}) where {N}
    return EmptyIterator{Vector{N}}()
end
