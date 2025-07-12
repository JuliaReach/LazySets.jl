"""
# Extended help

    intersection(LS1::LineSegment, LS2::LineSegment)

### Output

A `Singleton`, a `LineSegment`, or an `EmptySet` depending on the result of the
intersection.

### Notes

- If the line segments cross, or are parallel and have a single point in common,
  that point is returned.

- If the line segments are parallel and have a line segment in common, that
  segment is returned.

- Otherwise, there is no intersection and the empty set is returned.
"""
@validate function intersection(LS1::LineSegment, LS2::LineSegment)
    require(@__MODULE__, :LazySets; fun_name="intersection")

    # cast each segment as a line
    L1 = Line2D(LS1.p, LS1.q)
    L2 = Line2D(LS2.p, LS2.q)

    # find intersection between the lines
    m = intersection(L1, L2)
    N = promote_type(eltype(LS1), eltype(LS2))
    if m == L1
        # line segments are on the same line

        # check that the line segments intersect
        @inbounds begin
            # box approximation
            l, r = extrema((LS1.p[1], LS1.q[1]))
            b, t = extrema((LS1.p[2], LS1.q[2]))
            # check that the other line segment has at least one point inside the box
            if !(l <= LS2.p[1] <= r && b <= LS2.p[2] <= t) &&
               !(l <= LS2.q[1] <= r && b <= LS2.q[2] <= t)
                return EmptySet{N}(2)
            end
        end

        # find the middle two points of the four end points
        points = [LS1.p, LS1.q, LS2.p, LS2.q]
        sorted_points = sort(points; by=p -> (p[1], p[2]))
        @inbounds begin
            mid_point1 = sorted_points[2]
            mid_point2 = sorted_points[3]
        end
        if _isapprox(mid_point1, mid_point2)
            return Singleton(mid_point1)  # only one point in common
        else
            return LineSegment(mid_point1, mid_point2)
        end
    elseif m isa Singleton && m.element ∈ LS1 && m.element ∈ LS2
        return m  # the intersection point between the lines is in the segments
    else
        return EmptySet{N}(2)  # no intersection
    end
end
