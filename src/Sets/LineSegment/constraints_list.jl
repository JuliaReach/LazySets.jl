"""
# Extended help

    constraints_list(L::LineSegment)

### Input

- `sort_constraints` -- (optional; default: `false`) if `true`, the constraints
                        will be sorted according to `⪯`

### Algorithm

``L`` is defined by 4 constraints. In this algorithm, the first and the third
constraints are obtained by `halfspace_right` and `halfspace_left`,
respectively, and the other two are obtained by considering a vector parallel to
the line segment passing through one of the vertices such that the constraints
are sorted in counter-clockwise order; however, they are not necessarily sorted
according to `⪯`. (Use the `sort_constraints` option to achieve that.)
"""
function constraints_list(L::LineSegment; sort_constraints::Bool=false)
    p, q = L.p, L.q
    d = @inbounds [p[2] - q[2], q[1] - p[1]]
    clist = [halfspace_right(L), halfspace_right(p, p + d),
             halfspace_left(L), halfspace_left(q, q + d)]
    if sort_constraints
        # result is not sorted according to `⪯` yet, so sort it
        clist = _sort_constraints(clist)
    end
    return clist
end
