"""
    center(X::LazySet, i::Int)

Compute the center of a centrally symmetric set in a given dimension.

### Input

- `X` -- centrally symmetric set
- `i` -- dimension

### Output

A real number representing the center of the set in the given dimension.
"""
function center(::LazySet, ::Int) end

"""
    center(X::LazySet)

Compute the center of a centrally symmetric set.

### Input

- `X` -- centrally symmetric set

### Output

A vector with the center, or midpoint, of `X`.
"""
function center(::LazySet) end
