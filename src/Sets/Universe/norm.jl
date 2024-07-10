"""
    norm(U::Universe, [p]::Real=Inf)

Return the norm of a universe.
It is the norm of the enclosing ball (of the given ``p``-norm) of minimal volume
that is centered in the origin.

### Input

- `U` -- universe
- `p` -- (optional, default: `Inf`) norm

### Output

An error.
"""
function norm(U::Universe, p::Real=Inf)
    return error("a universe does not have a norm")
end
