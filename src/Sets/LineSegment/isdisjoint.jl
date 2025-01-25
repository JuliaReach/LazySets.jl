"""
# Extended help

    isdisjoint(L1::LineSegment, L2::LineSegment, [witness]::Bool=false)

### Algorithm

The algorithm is inspired from [here](https://stackoverflow.com/a/565282), which
itself is the special 2D case of a 3D algorithm from [Goldman90](@citet).

We first check if the two line segments are parallel, and if so, if they are
collinear. In the latter case, we check membership of any of the end points in
the other line segment. Otherwise the lines are not parallel, so we can solve an
equation of the intersection point, if it exists.
"""
function isdisjoint(L1::LineSegment, L2::LineSegment, witness::Bool=false)
    r = L1.q - L1.p
    if all(isapproxzero, r)
        # first line segment is a point
        empty_intersection = L1.q ∉ L2
        return _witness_result_empty(witness, empty_intersection, L1, L2, L1.q)
    end

    s = L2.q - L2.p
    if all(isapproxzero, s)
        # second line segment is a point
        empty_intersection = L2.q ∉ L1
        return _witness_result_empty(witness, empty_intersection, L1, L2, L2.q)
    end

    p1p2 = L2.p - L1.p
    u_numerator = right_turn(p1p2, r)
    u_denominator = right_turn(r, s)

    if u_denominator == 0
        # line segments are parallel
        if u_numerator == 0
            # line segments are collinear
            if L1.p ∈ L2
                empty_intersection = false
                if witness
                    v = L1.p
                end
            elseif L1.q ∈ L2
                empty_intersection = false
                if witness
                    v = L1.q
                end
            else
                empty_intersection = true
            end
        else
            # line segments are parallel and not collinear
            empty_intersection = true
        end
    else
        # line segments are not parallel
        u = u_numerator / u_denominator
        if u < 0 || u > 1
            empty_intersection = true
        else
            t = right_turn(p1p2, s) / u_denominator
            empty_intersection = t < 0 || t > 1
            if witness
                v = L1.p + t * r
            end
        end
    end
    if witness && !empty_intersection
        return (false, v)
    end
    return _witness_result_empty(witness, empty_intersection, L1, L2)
end
