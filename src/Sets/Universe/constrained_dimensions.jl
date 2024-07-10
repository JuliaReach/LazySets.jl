"""
    constrained_dimensions(U::Universe)

Return the indices in which a universe is constrained.

### Input

- `U` -- universe

### Output

The empty vector, as the universe is unconstrained in every dimension.
"""
function constrained_dimensions(U::Universe)
    return Int[]
end
