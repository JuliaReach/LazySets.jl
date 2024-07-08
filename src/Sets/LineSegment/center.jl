"""
    center(L::LineSegment)

Return the center of a 2D line segment.

### Input

- `L` -- 2D line segment

### Output

The center of the line segment.
"""
function center(L::LineSegment)
    return L.p + (L.q - L.p) / 2
end
