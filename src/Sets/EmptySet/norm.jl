"""
    norm(∅::EmptySet, [p]::Real=Inf)

Return the norm of an empty set.
It is the norm of the enclosing ball (of the given ``p``-norm) of minimal volume
that is centered in the origin.

### Input

- `∅` -- empty set
- `p` -- (optional, default: `Inf`) norm

### Output

An error.
"""
function norm(::EmptySet, ::Real=Inf)
    return error("the norm of an empty set is undefined")
end
