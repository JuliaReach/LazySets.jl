"""
    reflect(B::BallInf)

Concrete reflection of a ball in the infinity norm `B`, resulting in the
reflected set `-B`.

### Input

- `B` -- ball in the infinity norm

### Output

The `BallInf` representing `-B`.

### Algorithm

If ``B`` has center ``c`` and radius ``r``, then ``-B`` has center ``-c`` and
radius ``r``.
"""
function reflect(B::BallInf)
    return BallInf(-center(B), B.radius)
end
