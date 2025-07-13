"""
# Extended help

    σ(d::AbstractVector, L::LineSegment)

### Algorithm

If the angle between the vector ``q - p`` and ``d`` is bigger than 90° and less
than 270° (measured in counter-clockwise order), the result is ``p``, otherwise
it is ``q``.
If the angle is exactly 90° or 270°, or if the direction has norm zero, this
implementation returns ``q``.
"""
@validate function σ(d::AbstractVector, L::LineSegment)
    return sign(dot(L.q - L.p, d)) >= 0 ? L.q : L.p
end
