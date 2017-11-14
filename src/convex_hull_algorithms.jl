"""
    convex_hull(points; algorithm)

Compute the convex hull of points in the plane.

### Input

- `points`    -- array of vectors containing the 2D coordinates of the points
- `algorithm` -- (optional, default: `"monotone_chain"`) choose the convex
                 hull algorithm, valid options are:

    * `"monotone_chain"`
"""
function convex_hull(points; algorithm="monotone_chain")
    convex_hull!(copy(points), algorithm=algorithm)
end

"""
    convex_hull!(points; algorithm)

Compute the convex hull of points in the plane, in-place.
See also: `convex_hull`.
"""
function convex_hull!(points; algorithm="monotone_chain")
    if algorithm == "monotone_chain"
        return monotone_chain!(points)
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
    monotone_chain!(points)

Compute the convex hull of points in the plane using Andrew's monotone chain method.

### Input

- `points` -- array of vectors containing the 2D coordinates of the points;
              is sorted in-place inside this function

### Output

Array of vectors containing the 2D coordinates of the corner points of the
convex hull.

### Algorithm

This function implements Andrew's monotone chain convex hull algorithm
to construct the convex hull of a set of ``n`` points in the plane
in ``O(n \\log n)`` time.
For further details see the wikipedia page:
[Monotone chain](https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain)
"""
function monotone_chain!(points::Vector{S}) where{S<:AbstractVector{T}} where{T<:Real}

    @inline function build_hull!(semihull, iterator, points, zero_T)
        @inbounds for i in iterator
            while length(semihull) >= 2 && right_turn(semihull[end-1], semihull[end], points[i]) <= zero_T
                pop!(semihull)
            end
            push!(semihull, points[i])
        end
    end

    # sort the rows lexicographically (which requires a two-dimensional array)
    # points = sortrows(hcat(points...)', alg=QuickSort)  # out-of-place version
    sort!(points, by=x->(x[1], x[2]))                     # in-place version

    zero_T = zero(T)

    # build lower hull
    lower = Vector{S}()
    build_hull!(lower, indices(points)[1], points, zero_T)

    # build upper hull
    upper = Vector{S}()
    build_hull!(upper, reverse(indices(points)[1]), points, zero_T)

    # remove the last point of each segment because they are repeated
    return [lower[1:end-1]; upper[1:end-1]]
end
