"""
    norm(X::LazySet, [p]::Real=Inf)

Return the norm of a set.

### Input

- `X` -- set
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the norm.

### Notes

The norm of a set is the norm of the enclosing ball (of the given ``p``-norm) of
minimal volume that is centered in the origin.
"""
function norm(::LazySet, ::Real=Inf) end
