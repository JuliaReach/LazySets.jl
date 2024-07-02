"""
    reflect(B::Ball1)

Concrete reflection of a ball in the 1-norm `B`, resulting in the reflected set
`-B`.

### Input

- `B` -- ball in the 1-norm

### Output

The `Ball1` representing `-B`.

### Algorithm

If ``B`` has center ``c`` and radius ``r``, then ``-B`` has center ``-c`` and
radius ``r``.
"""
function reflect(B::Ball1)
    return Ball1(-center(B), B.radius)
end
