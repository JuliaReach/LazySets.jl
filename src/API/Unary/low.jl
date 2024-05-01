"""
    low(X::LazySet, i::Int)

Compute the lowest coordinate of a set in a given dimension.

### Input

- `X` -- set
- `i` -- dimension

### Output

A real number representing the lowest coordinate of the set in the given
dimension.

### Notes

The resulting value is the lower end of the projection.
"""
function low(::LazySet, ::Int) end

"""
    low(X::LazySet)

Compute the lowest coordinates of a set in each dimension.

### Input

- `X` -- set

### Output

A vector with the lowest coordinate of the set in each dimension.

### Notes

See also `low(X::LazySet, i::Int)`.

The result is the lowermost corner of the box approximation, so it is not
necessarily contained in `X`.
"""
function low(::LazySet) end
