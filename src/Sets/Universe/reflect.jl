"""
    reflect(U::Universe)

Concrete reflection of a universe `U`, resulting in the reflected set `-U`.

### Input

- `U` -- universe

### Output

The same universe.
"""
function reflect(U::Universe)
    return U
end
