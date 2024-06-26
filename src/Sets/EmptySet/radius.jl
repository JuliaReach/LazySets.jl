"""
    radius(∅::EmptySet, [p]::Real=Inf)

Return the radius of an empty set.
It is the radius of the enclosing ball (of the given ``p``-norm) of minimal
volume with the same center.

### Input

- `∅` -- empty set
- `p` -- (optional, default: `Inf`) norm

### Output

An error.
"""
function radius(::EmptySet, ::Real=Inf)
    return error("the radius of an empty set is undefined")
end
