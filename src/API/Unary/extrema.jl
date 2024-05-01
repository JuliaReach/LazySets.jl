"""
    extrema(X::LazySet, i::Int)

Compute the lowest and highest coordinate of a set in a given dimension.

### Input

- `X` -- set
- `i` -- dimension

### Output

Two real numbers representing the lowest and highest coordinate of the set in
the given dimension.

### Notes

The result is equivalent to `(low(X, i), high(X, i))`, but sometimes it can be
computed more efficiently.

The resulting values are the lower and upper ends of the projection.
"""
function extrema(::LazySet, ::Int) end

"""
    extrema(X::LazySet)

Compute the lowest and highest coordinate of a set in each dimension.

### Input

- `X` -- set

### Output

Two vectors with the lowest and highest coordinates of `X` in each dimension.

### Notes

See also [`extrema(X::LazySet, i::Int)`](@ref).

The result is equivalent to `(low(X), high(X))`, but sometimes it can be
computed more efficiently.

The resulting points are the lowest and highest corners of the box
approximation, so they are not necessarily contained in `X`.

### Algorithm

The default implementation computes the extrema via `low` and `high`.
"""
function extrema(::LazySet) end
