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

"""
    right_turn(O, A, B)

Determine if the acute angle defined by the three points O, A, B in the plane is
a right turn (counter-clockwise) with respect to the center O.

### Input

- `O` -- center point
- `A` -- one point
- `B` -- another point

### Algorithm

The [cross product](https://en.wikipedia.org/wiki/Cross_product) is used to
determine the sense of rotation. If the result is 0, the points are collinear;
if it is positive, the three points constitute a positive angle of rotation
around O from A to B; otherwise a negative angle.
"""
@inline right_turn(O, A, B) = (A[1] - O[1])*(B[2]-O[2]) - (A[2] - O[2])*(B[1]-O[1])

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

    # sort the rows lexicographically (which requires a two-dimensional array)
    points = sortrows(hcat(points...)')

    # build lower hull
    lower = Vector{Vector{T}}()
    for i in indices(points, 1)
        p = view(points, i, :)
        while length(lower) >= 2 && right_turn(lower[end-1], lower[end], p) <= zero(T)
            pop!(lower)
        end
        push!(lower, p)
    end

    # build upper hull
    upper = Vector{Vector{T}}()
    for i in size(points, 1):-1:1
        p = view(points, i, :)
        while length(upper) >= 2 && right_turn(upper[end-1], upper[end], p) <= zero(T)
            pop!(upper)
        end
        push!(upper, p)
    end

    # remove the last point of each segment because they are repeated
    return [lower[1:end-1]; upper[1:end-1]]
end
