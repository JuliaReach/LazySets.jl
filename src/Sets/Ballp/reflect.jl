"""
    reflect(B::Ballp)

Concrete reflection of a ball in the p-norm `B`, resulting in the reflected set
`-B`.

### Input

- `B` -- ball in the p-norm

### Output

The `Ballp` representing `-B`.

### Algorithm

If ``B`` has center ``c`` and radius ``r``, then ``-B`` has center ``-c`` and
radius ``r``. The norm remains the same.
"""
function reflect(B::Ballp)
    return Ballp(B.p, -center(B), B.radius)
end
