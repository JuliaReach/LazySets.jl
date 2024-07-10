"""
    radius(U::Universe, [p]::Real=Inf)

Return the radius of a universe.
It is the radius of the enclosing ball (of the given ``p``-norm) of minimal
volume.

### Input

- `U` -- universe
- `p` -- (optional, default: `Inf`) norm

### Output

An error.
"""
function radius(U::Universe, p::Real=Inf)
    return error("a universe does not have a radius")
end
