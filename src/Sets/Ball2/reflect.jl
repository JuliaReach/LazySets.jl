"""
    reflect(B::Ball2)

Concrete reflection of a ball in the 2-norm `B`, resulting in the reflected set
`-B`.

### Input

- `B` -- ball in the 2-norm

### Output

The `Ball2` representing `-B`.

### Algorithm

If ``B`` has center ``c`` and radius ``r``, then ``-B`` has center ``-c`` and
radius ``r``.
"""
function reflect(B::Ball2)
    return Ball2(-center(B), B.radius)
end
