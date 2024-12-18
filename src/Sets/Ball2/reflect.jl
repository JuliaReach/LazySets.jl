"""
# Extended help

    reflect(B::Ball2)

### Algorithm

If ``B`` has center ``c`` and radius ``r``, then ``-B`` has center ``-c`` and
radius ``r``.
"""
function reflect(B::Ball2)
    return Ball2(-center(B), B.radius)
end
