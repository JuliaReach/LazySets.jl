# Algorithm:
# Choose an arbitrary ray through `x`` and count whether the number of
# intersections with edges is odd. We choose the ray that is vertical upward.
# Vertical line segments are ignored (after checking whether `x` is a member).
@validate function in(x::AbstractVector, P::Polygon)
    require(@__MODULE__, :LazySets; fun_name="in")

    vlist = P.vertices
    if length(vlist) < 2
        if isempty(vlist)
            return false
        elseif length(vlist) == 1
            return @inbounds x == vlist[1]
        end
    end

    # vertical ray as a vertical line (compare y coordinates later)
    N = promote_type(eltype(x), eltype(P))
    vline = Line2D([one(N), zero(N)], @inbounds x[1])

    p = @inbounds vlist[end]
    odd = false
    @inbounds for q in vlist
        if (_leq(p[1], x[1]) && _geq(q[1], x[1])) || (_leq(q[1], x[1]) && _geq(p[1], x[1]))
            # line segment pq intersects vertical line through x
            if _isapprox(p[1], q[1])
                # vertical line segment
                if (_leq(p[2], x[2]) && _geq(q[2], x[2])) || (_leq(q[2], x[2]) && _geq(p[2], x[2]))
                    # x is on the line segment
                    return true
                end
            else
                # non-vertical line segment -> intersect (Line2D intersection is used)
                line2 = Line2D(p, q)
                y = element(intersection(line2, vline))
                # compare y coordinate
                if _geq(y[2], x[2])
                    if y == x
                        # x is on line segment
                        return true
                    end
                    odd = !odd
                end
            end
        end

        p = q
    end
    return odd
end
