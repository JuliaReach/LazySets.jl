"""
    andrew_monotone_chain(points)

### Input

- `points` -- array of vectors containing the 2D coordinates of the points

### Algorithm

This function implements Andrew's monotone chain convex hull algorithm, that constructs
the convex hull of a set of ``n`` 2D points in ``O(n log n)`` time.
For further details see:
[Monotone chain](https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain)
"""
function andrew_monotone_chain(points::Vector{Vector{T}}) where {T<:AbstractFloat}

    points = [sortrows(hcat(points...)')[i, :] for i in 1:length(points)]

    # Build lower hull
    lower = Vector{Vector{Float64}}()
    for p in points
        while length(lower) >= 2 && cross([lower[end]-lower[end-1]; 0], [p-lower[end-1]; 0])[3] <= 0
            pop!(lower)
        end
        push!(lower, p)
    end

    # Build upper hull
    upper = Vector{Vector{Float64}}()
    for i in length(points):-1:1
        p = points[i]
        while length(upper) >= 2 && cross([upper[end]-upper[end-1]; 0], [p-upper[end-1]; 0])[3] <= 0
            pop!(upper)
        end
        push!(upper, p)
    end

    return [lower[1:end-1]; upper[1:end-1]]
end
