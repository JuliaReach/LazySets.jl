"""
# Extended help

    an_element(L::LineSegment)

### Algorithm

The output is the first vertex of the line segment.
"""
@validate function an_element(L::LineSegment)
    return L.p
end
