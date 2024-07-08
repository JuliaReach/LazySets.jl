"""
    dim(L::LineSegment)

Return the ambient dimension of a 2D line segment.

### Input

- `L` -- 2D line segment

### Output

The ambient dimension of the 2D line segment, which is ``2``.
"""
function dim(L::LineSegment)
    return 2
end
