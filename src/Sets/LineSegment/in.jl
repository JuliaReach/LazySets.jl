"""
# Extended help

    ∈(x::AbstractVector, L::LineSegment)

### Algorithm

Let ``L = (p, q)`` be the line segment with extreme points ``p`` and ``q``, and
let ``x`` be the given point.

1. A necessary condition for ``x ∈ (p, q)`` is that the three points are
   aligned, thus their cross product should be zero.
2. It remains to check that ``x`` belongs to the box approximation of ``L``.
   This amounts to comparing each coordinate with those of the extremes ``p``
   and ``q``.

### Notes

The algorithm is inspired from [here](https://stackoverflow.com/a/328110).
"""
@validate function ∈(x::AbstractVector, L::LineSegment)
    # check if point x is on the line through the line segment (p, q)
    p = L.p
    q = L.q
    if isapproxzero(right_turn(p, q, x))
        # check if the point is inside the box approximation of the line segment
        return @inbounds (_leq(min(p[1], q[1]), x[1]) &&
                          _leq(x[1], max(p[1], q[1])) &&
                          _leq(min(p[2], q[2]), x[2]) &&
                          _leq(x[2], max(p[2], q[2])))
    else
        return false
    end
end
