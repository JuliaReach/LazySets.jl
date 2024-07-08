"""
    halfspace_left(L::LineSegment)

Return a half-space describing the 'left' of a two-dimensional 2D line segment
through two points.

### Input

- `L` -- 2D line segment

### Output

The half-space whose boundary goes through the two points `p` and `q` and which
describes the left-hand side of the directed line segment `pq`.
"""
halfspace_left(L::LineSegment) = halfspace_left(L.p, L.q)
