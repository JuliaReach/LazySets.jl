function default_convex_hull_algorithm(points)
    if length(points) != 0 && length(first(points)) == 2
        return "monotone_chain"
    else
        return nothing
    end
end

"""
    convex_hull(points::Vector{VN};
                [algorithm]=default_convex_hull_algorithm(points),
                [backend]=nothing
                )::Vector{VN} where {N<:Real, VN<:AbstractVector{N}}

Compute the convex hull of the given points.

### Input

- `points`    -- list of vectors
- `algorithm` -- (optional, default: depends on the dimension) the convex hull
                 algorithm, see valid options below
- `backend`   -- (optional, default: `"nothing"`) polyhedral computation backend
                 for higher-dimensional point sets

### Output

The convex hull as a list of vectors with the coordinates of the points.

### Algorithm

A pre-processing step treats the cases with `0`, `1` and `2` points in any dimension.
For more than `3` points, the algorithm used depends on the dimension.

For the one-dimensional case we return the minimum and maximum points, in that order.

The two-dimensional case is handled with a planar convex hull algorithm. The
following algorithms are available:

* `"monotone_chain"`        -- compute the convex hull of points in the plane
                               using Andrew's monotone chain method
* `"monotone_chain_sorted"` -- the same as `"monotone_chain"` but assuming that
                               the points are already sorted in counter-clockwise
                               fashion

See the reference docstring of each of those algorithms for details.

The higher dimensional case is treated using the concrete polyhedra library
`Polyhedra`, that gives access to libraries such as `CDDLib` and `ConvexHull.jl`.
These libraries can be chosen from the `backend` argument.

### Notes

For the in-place version use `convex_hull!` instead of `convex_hull`.

### Examples

Compute the convex hull of a random set of points:

```jldoctest ch_label
julia> points = [randn(2) for i in 1:30]; # 30 random points in 2D

julia> hull = convex_hull(points);

julia> typeof(hull)
Array{Array{Float64,1},1}
```

Plot both the random points and the computed convex hull polygon:

```jldoctest ch_label
julia> using Plots;

julia> plot([Tuple(pi) for pi in points], seriestype=:scatter);

julia> plot!(VPolygon(hull), alpha=0.2);
```
"""
function convex_hull(points::Vector{VN};
                     algorithm=default_convex_hull_algorithm(points),
                     backend=nothing
                     )::Vector{VN} where {N<:Real, VN<:AbstractVector{N}}
    return convex_hull!(copy(points), algorithm=algorithm, backend=backend)
end

function convex_hull!(points::Vector{VN};
                      algorithm=default_convex_hull_algorithm(points),
                      backend=nothing)::Vector{VN} where {N<:Real, VN<:AbstractVector{N}}

    m = length(points)

    # =====================
    # Treat special cases
    # =====================

    # zero or one point
    if m == 0 || m == 1
        return points
    end

    n = length(first(points))

    # two points
    if m == 2
        if n == 1
            p1, p2 = points[1], points[2]
            if p1 == p2  # check for redundancy
                pop!(points)
            elseif p1[1] > p2[1]
                points[1], points[2] = p2, p1
            end
        elseif n == 2
            # special case, see #876
            p1, p2 = points[1], points[2]
            if p1 == p2  # check for redundancy
                pop!(points)
            elseif p1 > p2
                points[1], points[2] = p2, p1
            end
        end
        return points
    end

    # ===========================
    # Dispatch for general cases
    # ===========================

    # _convex_hull_x can assume that there are at least three points
    if n == 1
        return _convex_hull_1d!(points)
    elseif n == 2 
        return _convex_hull_2d!(points, algorithm=algorithm)
    else
        return _convex_hull_nd!(points, backend=backend)
    end
end

function _convex_hull_1d!(points::Vector{VN})::Vector{VN} where {N<:Real, VN<:AbstractVector{N}}
    points[1:2] = [minimum(points), maximum(points)]
    return resize!(points, 2)
end

function _convex_hull_nd!(points::Vector{VN};
                          backend=nothing
                          )::Vector{VN} where {N<:Real, VN<:AbstractVector{N}}
    V = VPolytope(points)
    if backend == nothing
        backend = default_polyhedra_backend(V, N)
    end
    Vch = remove_redundant_vertices(V, backend=backend)
    m = length(Vch.vertices)
    points[1:m] = Vch.vertices
    return resize!(points, m)
end

function _convex_hull_2d!(points::Vector{VN};
                          algorithm::String="monotone_chain"
                          )::Vector{VN} where {N<:Real, VN<:AbstractVector{N}}
    if algorithm == "monotone_chain"
        return monotone_chain!(points)
    elseif algorithm == "monotone_chain_sorted"
        return monotone_chain!(points, sort=false)
    else
        error("the convex hull algorithm $algorithm is unknown")
    end
end

"""
    right_turn(O::AbstractVector{N}, A::AbstractVector{N}, B::AbstractVector{N}
              )::N where {N<:Real}

Determine if the acute angle defined by the three points `O`, `A`, `B` in the
plane is a right turn (counter-clockwise) with respect to the center `O`.

### Input

- `O` -- 2D center point
- `A` -- 2D one point
- `B` -- 2D another point

### Output

Scalar representing the rotation.

### Algorithm

The [cross product](https://en.wikipedia.org/wiki/Cross_product) is used to
determine the sense of rotation. If the result is 0, the points are collinear;
if it is positive, the three points constitute a positive angle of rotation
around `O` from `A` to `B`; otherwise they constitute a negative angle.
"""
@inline function right_turn(O::AbstractVector{N},
                            A::AbstractVector{N},
                            B::AbstractVector{N})::N where {N<:Real}
    return (A[1] - O[1]) * (B[2] - O[2]) - (A[2] - O[2]) * (B[1] - O[1])
end

"""
    monotone_chain!(points::Vector{VN}; sort::Bool=true
                   )::Vector{VN} where {N<:Real, VN<:AbstractVector{N}}

Compute the convex hull of points in the plane using Andrew's monotone chain
method.

### Input

- `points` -- list of 2D vectors; is sorted in-place inside this function
- `sort`   -- (optional, default: `true`) flag for sorting the vertices
              lexicographically; sortedness is required for correctness

### Output

List of vectors containing the 2D coordinates of the corner points of the
convex hull.

### Notes

For large sets of points, it is convenient to use static vectors to get
maximum performance. For information on how to convert usual vectors
into static vectors, see the type `SVector` provided by the
[StaticArrays](http://juliaarrays.github.io/StaticArrays.jl/stable/)
package.

### Algorithm

This function implements Andrew's monotone chain convex hull algorithm to
construct the convex hull of a set of ``n`` points in the plane in
``O(n \\log n)`` time.
For further details see
[Monotone chain](https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain)
"""
function monotone_chain!(points::Vector{VN}; sort::Bool=true
                        )::Vector{VN} where {N<:Real, VN<:AbstractVector{N}}

    @inline function build_hull!(semihull, iterator, points, zero_N)
        @inbounds for i in iterator
            while length(semihull) >= 2 &&
                    (right_turn(semihull[end-1], semihull[end], points[i])
                         <= zero_N)
                pop!(semihull)
            end
            push!(semihull, points[i])
        end
    end

    if sort
        # sort the rows lexicographically (requires a two-dimensional array)
        # points = sortrows(hcat(points...)', alg=QuickSort) # out-of-place version
        sort!(points, by=x->(x[1], x[2]))                    # in-place version
    end

    zero_N = zero(N)

    # build lower hull
    lower = Vector{VN}()
    build_hull!(lower, axes(points)[1], points, zero_N)

    # build upper hull
    upper = Vector{VN}()
    build_hull!(upper, reverse(axes(points)[1]), points, zero_N)

    # remove the last point of each segment because they are repeated
    copyto!(points, @view(lower[1:end-1]))
    copyto!(points, length(lower), @view(upper[1:end-1]))
    return resize!(points, length(lower) + length(upper) - 2)
end
