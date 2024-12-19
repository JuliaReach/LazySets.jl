"""
# Extended help

    constraints_list(L::LineSegment)

### Algorithm

``L`` is defined by 4 constraints. In this algorithm, the first two constraints
are returned by ``halfspace_right`` and ``halfspace_left``, and the other two
are obtained by considering a vector parallel to the line segment passing
through one of the vertices.
"""
function constraints_list(L::LineSegment)
    p, q = L.p, L.q
    d = @inbounds [p[2] - q[2], q[1] - p[1]]
    return [halfspace_left(L), halfspace_right(L),
            halfspace_right(p, p + d), halfspace_left(q, q + d)]
end
