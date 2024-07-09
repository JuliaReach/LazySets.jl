"""
    vertices_list(L::LineSegment)

Return the list of vertices of a 2D line segment.

### Input

- `L` -- 2D line segment

### Output

The list of end points of the line segment.
"""
function vertices_list(L::LineSegment)
    return [L.p, L.q]
end
