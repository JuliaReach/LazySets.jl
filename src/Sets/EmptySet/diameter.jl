"""
    diameter(∅::EmptySet, [p]::Real=Inf)

Return the diameter of an empty set.
It is the maximum distance between any two elements of the set or, equivalently,
the diameter of the enclosing ball (of the given ``p``-norm) of minimal volume
with the same center.

### Input

- `∅` -- empty set
- `p` -- (optional, default: `Inf`) norm

### Output

An error.
"""
function diameter(::EmptySet, ::Real=Inf)
    return error("the diameter of an empty set is undefined")
end
