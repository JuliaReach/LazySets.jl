"""
    halfspace_right(L::LineSegment)

Return a half-space describing the 'right' of a two-dimensional 2D line segment
through two points.

### Input

- `L` -- 2D line segment

### Output

The half-space whose boundary goes through the two points `p` and `q` and which
describes the right-hand side of the directed line segment `pq`.
"""
halfspace_right(L::LineSegment) = halfspace_right(L.p, L.q)
