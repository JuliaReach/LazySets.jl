"""
    intersection(LS1::LineSegment, LS2::LineSegment)

Compute the intersection of two line segments.

### Input

- `LS1` -- line segment
- `LS2` -- line segment

### Output

A singleton, line segment, or the empty set depending on the result of the
intersection.

### Notes

- If the line segments cross, or are parallel and have one point in common,
  that point is returned.

- If the line segments are parallel and have a line segment in common, that
  segment is returned.

- Otherwise, if there is no intersection, an empty set is returned.
"""
function intersection(LS1::LineSegment, LS2::LineSegment)

    # cast each segment as a line
    L1 = Line2D(LS1.p, LS1.q)
    L2 = Line2D(LS2.p, LS2.q)

    # find intersection between the lines
    m = intersection(L1, L2)
    N = promote_type(eltype(LS1), eltype(LS2))
    if m == L1
        # determine which segment is in both
        p1 = max(min(LS1.p[1], LS1.q[1]), min(LS2.p[1], LS2.q[1]))
        p2 = max(min(LS1.p[2], LS1.q[2]), min(LS2.p[2], LS2.q[2]))
        q1 = min(max(LS1.p[1], LS1.q[1]), max(LS2.p[1], LS2.q[1]))
        q2 = min(max(LS1.p[2], LS1.q[2]), max(LS2.p[2], LS2.q[2]))
        if _isapprox(p1, q1) && _isapprox(p2, q2)
            return Singleton([p1, p2])  # edges have a point in common

        elseif _leq(p1, q1) && _leq(p2, q2)
            return LineSegment([p1, p2], [q1, q2])

        else
            return EmptySet{N}(2)  # no intersection
        end

    elseif m isa Singleton && m.element ∈ LS1 && m.element ∈ LS2
        return m  # the intersection point between the lines is in the segments

    else
        return EmptySet{N}(2)  # no intersection
    end
end
