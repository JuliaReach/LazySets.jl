"""
    convex_hull(points; algorithm)

Compute the convex hull of points in the plane.

### Input

- `points` -- array of vectors containing the 2D coordinates of the points

"""
function convex_hull(points; algorithm="andrew")
    if algorithm == "andrew"
        return andrew_monotone_chain(points)
    else
        error("this convex hull algorithm is unknown")
    end
end

#=
@inline function right_turn(O, A, B)

end
=#

"""
    andrew_monotone_chain(points)

Compute the convex hull of points in the plane using Andrew's monotone chain method.

### Input

- `points` -- array of vectors containing the 2D coordinates of the points

### Algorithm

This function implements Andrew's monotone chain convex hull algorithm
to construct the convex hull of a set of ``n`` points in the plane
in ``O(n \\log n)`` time.
For further details see the wikipedia page:
[Monotone chain](https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain)
"""
function andrew_monotone_chain(points::Vector{Vector{T}}) where {T<:Real}

    points = [sortrows(hcat(points...)')[i, :] for i in 1:length(points)]

    # Build lower hull
    lower = Vector{Vector{T}}()
    for p in points
        while length(lower) >= 2 && cross([lower[end]-lower[end-1]; zero(T)], [p-lower[end-1]; zero(T)])[3] <= zero(T)
            pop!(lower)
        end
        push!(lower, p)
    end

    # Build upper hull
    upper = Vector{Vector{T}}()
    for i in length(points):-1:1
        p = points[i]
        while length(upper) >= 2 && cross([upper[end]-upper[end-1]; zero(T)], [p-upper[end-1]; zero(T)])[3] <= zero(T)
            pop!(upper)
        end
        push!(upper, p)
    end

    return [lower[1:end-1]; upper[1:end-1]]
end
