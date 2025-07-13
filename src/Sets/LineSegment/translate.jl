"""
# Extended help

    translate(L::LineSegment, v::AbstractVector)

### Algorithm

We add the vector to both defining points of the line segment.
"""
@validate function translate(L::LineSegment, v::AbstractVector)
    return LineSegment(L.p + v, L.q + v)
end
