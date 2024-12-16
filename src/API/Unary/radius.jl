"""
    radius(X::LazySet, [p]::Real=Inf)

Return the radius of a set.

### Input

- `X` -- set
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the radius.

### Notes

The radius of a set is the radius of the enclosing ball (of the given
``p``-norm) of minimal volume.
"""
function radius(::LazySet, ::Real=Inf) end
