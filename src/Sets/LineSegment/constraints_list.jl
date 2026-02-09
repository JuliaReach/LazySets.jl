"""
# Extended help

    constraints_list(L::LineSegment)

### Algorithm

``L`` is defined by 4 constraints. In this algorithm, the first and the third
constraints are obtained by `halfspace_right` and `halfspace_left`,
respectively, and the other two are obtained by considering a vector parallel to
the line segment passing through one of the vertices such that the constraints
are sorted in counter-clockwise order.
"""
function constraints_list(L::LineSegment)
    p, q = L.p, L.q
    d = @inbounds [p[2] - q[2], q[1] - p[1]]
    return [halfspace_right(L), halfspace_right(p, p + d),
            halfspace_left(L), halfspace_left(q, q + d)]
end
