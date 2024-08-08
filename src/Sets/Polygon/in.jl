# Algorithm:
# Choose an arbitrary ray through `x`` and count whether the number of
# intersections with edges is odd. We choose the ray that is vertical upward.
# Vertical line segments are ignored (after checking whether `x` is a member).
function ∈(x::AbstractVector, P::Polygon)
    require(@__MODULE__, :LazySets; fun_name="in")
    @assert length(x) == dim(P) "a vector of length $(length(x)) is " *
                                "incompatible with a $(dim(P))-dimensional set"

    N = promote_type(eltype(x), eltype(P))

    vlist = P.vertices
    if length(vlist) < 2
        if isempty(vlist)
            return false
        elseif length(vlist) == 1
            return @inbounds x == vlist[1]
        end
    end

    # vertical ray as a vertical line (compare y coordinates later)
    vline = Line2D([one(N), zero(N)], @inbounds x[1])

    p = @inbounds vlist[end]
    odd = false
    @inbounds for q in vlist
        if (p[1] <= x[1] && q[1] >= x[1]) || (q[1] <= x[1] && p[1] >= x[1])
            # line segment pq intersects vertical line through x
            if p[1] == q[1]
                # vertical line segment
                if (p[2] <= x[2] && q[2] >= x[2]) || (q[2] <= x[2] && p[2] >= x[2])
                    # x is on the line segment
                    return true
                end
            else
                # non-vertical line segment -> intersect (Line2D intersection is used)
                line2 = Line2D(p, q)
                y = element(intersection(line2, vline))
                # compare y coordinate
                if y[2] >= x[2]
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
