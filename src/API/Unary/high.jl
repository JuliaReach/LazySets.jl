"""
    high(X::LazySet, i::Int)

Compute the highest coordinate of a set in a given dimension.

### Input

- `X` -- set
- `i` -- dimension

### Output

A real number representing the highest coordinate of the set in the given
dimension.

### Notes

The resulting value is the upper end of the projection.
"""
function high(::LazySet, ::Int) end

"""
    high(X::LazySet)

Compute the highest coordinate of a set in each dimension.

### Input

- `X` -- set

### Output

A vector with the highest coordinate of the set in each dimension.

### Notes

See also `high(X::LazySet, i::Int)`.

The result is the uppermost corner of the box approximation, so it is not
necessarily contained in `X`.
"""
function high(::LazySet) end
