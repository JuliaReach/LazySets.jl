"""
    diameter(U::Universe, [p]::Real=Inf)

Return the diameter of a universe.
It is the diameter of the enclosing ball (of the given ``p``-norm) of minimal
volume .

### Input

- `U` -- universe
- `p` -- (optional, default: `Inf`) norm

### Output

An error.
"""
function diameter(U::Universe, p::Real=Inf)
    return error("a universe does not have a diameter")
end
