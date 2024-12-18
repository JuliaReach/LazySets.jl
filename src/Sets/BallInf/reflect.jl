"""
# Extended help

    reflect(B::BallInf)

### Algorithm

If ``B`` has center ``c`` and radius ``r``, then ``-B`` has center ``-c`` and
radius ``r``.
"""
function reflect(B::BallInf)
    return BallInf(-center(B), B.radius)
end
