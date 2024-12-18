"""
# Extended help

    reflect(B::Ball1)

### Algorithm

If ``B`` has center ``c`` and radius ``r``, then ``-B`` has center ``-c`` and
radius ``r``.
"""
function reflect(B::Ball1)
    return Ball1(-center(B), B.radius)
end
