"""
    diameter(X::LazySet, [p]::Real=Inf)

Return the diameter of a set.

### Input

- `X` -- set
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the diameter.

### Notes

The diameter of a set is the maximum distance between any two elements of the
set, or, equivalently, the diameter of the enclosing ball (of the given
``p``-norm) of minimal volume with the same center.
"""
function diameter(::LazySet, ::Real=Inf) end
