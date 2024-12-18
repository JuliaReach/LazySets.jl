"""
# Extended help

    reflect(B::Ballp)

### Algorithm

If ``B`` has center ``c`` and radius ``r``, then ``-B`` has center ``-c`` and
radius ``r``. The norm remains the same.
"""
function reflect(B::Ballp)
    return Ballp(B.p, -center(B), B.radius)
end
