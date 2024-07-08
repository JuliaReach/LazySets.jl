"""
    an_element(L::LineSegment)

Return some element of a 2D line segment.

### Input

- `L` -- 2D line segment

### Output

The first vertex of the line segment.
"""
function an_element(L::LineSegment)
    return L.p
end
