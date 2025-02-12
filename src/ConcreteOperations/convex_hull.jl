function default_convex_hull_algorithm(points)
    if length(points) != 0 && length(first(points)) == 2
        return "monotone_chain"
    else
        return nothing
    end
end

"""
    convex_hull(X::LazySet, Y::LazySet; [algorithm]=nothing,
                [backend]=nothing, [solver]=nothing)

Compute the convex hull of two polytopic sets.

### Input

- `X`         -- polytopic set
- `Y`         -- polytopic set
- `algorithm` -- (optional, default: `nothing`) the convex-hull algorithm
- `backend`   -- (optional, default: `nothing`) backend for polyhedral
                 computations (used for higher-dimensional sets)
- `solver`    -- (optional, default: `nothing`) the linear-programming solver
                 used in the backend

### Output

If the sets are empty, the result is an `EmptySet`.
If the convex hull consists of a single point, the result is a `Singleton`.
If the input sets are one-dimensional, the result is an `Interval`.
If the input sets are two-dimensional, the result is a `VPolygon`.
Otherwise the result is a `VPolytope`.

### Algorithm

We compute the vertices of both `X` and `Y` using `vertices_list` and then
compute the convex hull of the union of those vertices.
"""
function convex_hull(X::LazySet, Y::LazySet; algorithm=nothing,
                     backend=nothing, solver=nothing)
    n = dim(X)
    @assert n == dim(Y) "the convex hull requires two sets of the same " *
                        "dimension, but they have dimension $n and $(dim(Y))"

    vlist = convex_hull!([vertices_list(X); vertices_list(Y)];
                         algorithm=algorithm, backend=backend, solver=solver)
    return _convex_hull_set(vlist; n=n)
end

# turn a list of points describing the convex hull into an appropriate set
function _convex_hull_set(vlist; n::Int=(isempty(vlist) ? -1 : length(vlist[1])))
    m = length(vlist)
    if m == 0
        N = eltype(vlist)
        return EmptySet{N}(n)
    elseif m == 1
        return Singleton(vlist[1])
    elseif n == 1
        @assert m == 2 "a convex hull in 1D can have at most two vertices"
        # points are sorted
        return Interval(vlist[1][1], vlist[end][1])
    elseif n == 2
        # points are sorted
        return VPolygon(vlist)
    end
    return VPolytope(vlist)
end

"""
    convex_hull(points::Vector{VN}; [algorithm]=nothing, [backend]=nothing,
                [solver]=nothing) where {N, VN<:AbstractVector{N}}

Compute the convex hull of a list of points.

### Input

- `points`    -- list of points
- `algorithm` -- (optional, default: `nothing`) the convex-hull algorithm; see
                 below for valid options
- `backend`   -- (optional, default: `nothing`) polyhedral computation backend
                 for higher-dimensional point sets
- `solver`    -- (optional, default: `nothing`) the linear-programming solver
                 used in the backend

### Output

The convex hull as a list of points.

### Algorithm

A pre-processing step treats the cases with up to two points for one dimension
and up to four points for two dimensions.
For more points in one resp. two dimensions, we use more general algorithms.

For the one-dimensional case, we return the minimum and maximum points, in that
order.

The two-dimensional case is handled with a planar convex-hull algorithm. The
following algorithms are available:

* `"monotone_chain"`        -- compute the convex hull of points in the plane
                               using Andrew's monotone-chain method
* `"monotone_chain_sorted"` -- the same as `"monotone_chain"` but assuming that
                               the points are already sorted in
                               counter-clockwise fashion

See the reference docstring of each of those algorithms for details.

The higher-dimensional case is treated using the concrete polyhedra library
`Polyhedra`, which gives access to libraries such as `CDDLib` and
`ConvexHull.jl`. These libraries can be chosen via the `backend` argument.

### Notes

For the in-place version use `convex_hull!` instead of `convex_hull`.

### Examples

Compute the convex hull of a random set of points:

```jldoctest ch_label
julia> points = [randn(2) for i in 1:30]; # 30 random points in 2D

julia> hull = convex_hull(points);

julia> typeof(hull)
Vector{Vector{Float64}} (alias for Array{Array{Float64, 1}, 1})
```
"""
function convex_hull(points::Vector{VN}; algorithm=nothing, backend=nothing,
                     solver=nothing) where {N,VN<:AbstractVector{N}}
    return convex_hull!(copy(points); algorithm=algorithm, backend=backend,
                        solver=solver)
end

function convex_hull!(points::Vector{VN}; algorithm=nothing, backend=nothing,
                      solver=nothing) where {N,VN<:AbstractVector{N}}
    m = length(points)

    # zero or one point
    if m == 0 || m == 1
        return points
    end

    n = length(first(points))

    if n == 1 # dimensional check
        if m == 2
            # two points case in 1d
            return _two_points_1d!(points)
        else
            # general case in 1d
            return _convex_hull_1d!(points)
        end
    elseif n == 2
        return _convex_hull_2d_preprocess!(points, m; algorithm=algorithm)
    else
        # general case in nd
        return _convex_hull_nd!(points; backend=backend, solver=solver)
    end
end

function _convex_hull_2d_preprocess!(points, m=length(points); algorithm=nothing)
    if m == 1
        return points
    elseif m == 2
        # two points case in 2d
        return _two_points_2d!(points)
    elseif m == 3
        # three points case in 2d
        return _three_points_2d!(points)
    elseif m == 4
        # four points case in 2d
        return _four_points_2d!(points)
    else
        # general case in 2d
        return _convex_hull_2d!(points; algorithm=algorithm)
    end
end

function _two_points_1d!(points)
    p1, p2 = points[1], points[2]
    if _isapprox(p1[1], p2[1]) # check for redundancy
        pop!(points)
    elseif p1[1] > p2[1]
        points[1] = p2
        points[2] = p1
    end
    return points
end

function _two_points_2d!(points)
    p1, p2 = points[1], points[2]
    if _isapprox(p1[1], p2[1]) && _isapprox(p1[2], p2[2]) # check for redundancy
        pop!(points)
    end
    return points
end

function _three_points_2d!(points::AbstractVector{<:AbstractVector{N}}) where {N}
    # Algorithm: the function takes three points and uses the formula from here:
    # https://stackoverflow.com/a/2122620
    # to decide if the points are ordered in counter-clockwise order. The result
    # is saved in the 'turn' flag and returned in ccw fashion acting according
    # to 'turn'. For the cases where the points are collinear, we pass the
    # points with the minimum and maximum first component to the function for
    # two points in 2d(_two_points_2d). If those are equal, we do the same but
    # with the second component.
    A, B, C = points[1], points[2], points[3]
    turn = right_turn(A, B, C)

    if isapproxzero(turn)
        # ABC are collinear
        if isapprox(A[1], B[1]) && isapprox(B[1], C[1]) && isapprox(C[1], A[1])
            # points are approximately equal in their first component
            if isapprox(A[2], B[2]) && isapprox(B[2], C[2]) && isapprox(C[2], A[2])
                # all points are approximately equal
                pop!(points)
                pop!(points)
                return points
            else
                points[1] = points[argmin([A[2], B[2], C[2]])]
                points[2] = points[argmax([A[2], B[2], C[2]])]
            end
        else
            points[1] = points[argmin([A[1], B[1], C[1]])]
            points[2] = points[argmax([A[1], B[1], C[1]])]
        end
        pop!(points)
        _two_points_2d!(points)
    elseif turn < zero(N)
        # ABC is CW
        points[1], points[2], points[3] = C, B, A
    end
    # else ABC is CCW => nothing to do
    return points
end

# given an index i in {1, 2, 3}, return the element among (A, B, C) in the i-th
# position
@inline function _get_i(i, A, B, C)
    return i == 1 ? A : i == 2 ? B : C
end

function _collinear_case!(points, A, B, C, D)
    # A, B and C collinear, D is the extra point
    if isapprox(A[1], B[1]) && isapprox(B[1], C[1]) && isapprox(C[1], A[1])
        # points are approximately equal in their first component
        if isapprox(A[2], B[2]) && isapprox(B[2], C[2]) && isapprox(C[2], A[2])
            # the three points are approximately equal
            points[1], points[2] = A, D
            pop!(points)
            pop!(points)
            return _two_points_2d!(points)
        else
            # assign the points with max and min value in their second component
            # to the firsts points and the extra point to the third place, then
            # pop the point that was in the middle
            min_y, max_y = _arg_minmax(A[2], B[2], C[2])
            points[1] = _get_i(min_y, A, B, C)
            points[2] = _get_i(max_y, A, B, C)
            points[3] = D
            pop!(points)
        end
    else
        # assign the points with max and min value in their first component to
        # the firsts points and the extra point to the third place, then pop the
        # point that was in the middle
        min_x, max_x = _arg_minmax(A[1], B[1], C[1])
        points[1] = _get_i(min_x, A, B, C)
        points[2] = _get_i(max_x, A, B, C)
        points[3] = D
        pop!(points)
    end
    return _three_points_2d!(points)
end

# return the indices of the minimum and maximum of three numbers a, b, c
function _arg_minmax(a, b, c)
    if a > b
        min, max = b, a
        imin, imax = 2, 1
    else
        min, max = a, b
        imin, imax = 1, 2
    end
    if c > max
        imax = 3
    elseif c < min
        imin = 3
    end
    return imin, imax
end

function _four_points_2d!(points::AbstractVector{<:AbstractVector{N}}) where {N}
    A, B, C, D = points[1], points[2], points[3], points[4]
    tri_ABC = right_turn(A, B, C)
    tri_ABD = right_turn(A, B, D)
    tri_BCD = right_turn(B, C, D)
    tri_CAD = right_turn(C, A, D)
    key = 0
    if tri_ABC > zero(N)
        key = key + 1000
    end
    if tri_ABD > zero(N)
        key = key + 100
    end
    if tri_BCD > zero(N)
        key = key + 10
    end
    if tri_CAD > zero(N)
        key = key + 1
    end

    if isapproxzero(tri_ABC)
        return _collinear_case!(points, A, B, C, D)
    elseif isapproxzero(tri_ABD)
        return _collinear_case!(points, A, B, D, C)
    elseif isapproxzero(tri_BCD)
        return _collinear_case!(points, B, C, D, A)
    elseif isapproxzero(tri_CAD)
        return _collinear_case!(points, C, A, D, B)
    end

    # ABC  ABD  BCD  CAD  hull
    # ------------------------
    #  +    +    +    +   ABC
    #  +    +    +    -   ABCD
    #  +    +    -    +   ABDC
    #  +    +    -    -   ABD
    #  +    -    +    +   ADBC
    #  +    -    +    -   BCD
    #  +    -    -    +   CAD
    #  +    -    -    -   [should not happen]
    #  -    +    +    +   [should not happen]
    #  -    +    +    -   ACD
    #  -    +    -    +   DCB
    #  -    +    -    -   DACB
    #  -    -    +    +   ADB
    #  -    -    +    -   ACDB
    #  -    -    -    +   ADCB
    #  -    -    -    -   ACB
    if key == 1111
        points[1], points[2], points[3] = A, B, C #  +    +    +    +   ABC
        pop!(points)
    elseif key == 1110
        points[1], points[2], points[3], points[4] = A, B, C, D #  +    +    +    -   ABCD
    elseif key == 1101
        points[1], points[2], points[3], points[4] = A, B, D, C #  +    +    -    +   ABDC
    elseif key == 1100
        points[1], points[2], points[3] = A, B, D               #  +    +    -    -   ABD
        pop!(points)
    elseif key == 1011
        points[1], points[2], points[3], points[4] = A, D, B, C #  +    -    +    +   ADBC
    elseif key == 1010
        points[1], points[2], points[3] = B, C, D               #  +    -    +    -   BCD
        pop!(points)
    elseif key == 1001
        points[1], points[2], points[3] = C, A, D               #  +    -    -    +    CAD
        pop!(points)
    elseif key == 0110
        points[1], points[2], points[3] = A, C, D               #  -    +    +    -   ACD
        pop!(points)
    elseif key == 0101
        points[1], points[2], points[3] = D, C, B               #  -    +    -    +   DCB
        pop!(points)
    elseif key == 0100
        points[1], points[2], points[3], points[4] = D, A, C, B #  -    +    -    -   DACB
    elseif key == 0011
        points[1], points[2], points[3] = A, D, B               #  -    -    +    +   ADB
        pop!(points)
    elseif key == 0010
        points[1], points[2], points[3], points[4] = A, C, D, B #  -    -    +    -   ACDB
    elseif key == 0001
        points[1], points[2], points[3], points[4] = A, D, C, B #  -    -    -    +   ADCB
    elseif key == 0000
        points[1], points[2], points[3] = A, C, B               #  -    -    -    -   ACB
        pop!(points)
    else
        error("unexpected case in convex-hull algorithm")
    end
    return points
end

function _convex_hull_1d!(points::Vector{VN}) where {N,VN<:AbstractVector{N}}
    points[1:2] = [minimum(points), maximum(points)]
    return resize!(points, 2)
end

function _convex_hull_nd!(points::Vector{VN};
                          backend=nothing,
                          solver=nothing) where {N,VN<:AbstractVector{N}}
    V = VPolytope(points)
    Vch = VPolytopeModule.remove_redundant_vertices(V; backend=backend, solver=solver)
    m = length(Vch.vertices)
    points[1:m] = Vch.vertices
    return resize!(points, m)
end

function _convex_hull_2d!(points::Vector{VN};
                          algorithm="monotone_chain") where {N,VN<:AbstractVector{N}}
    if isnothing(algorithm)
        algorithm = default_convex_hull_algorithm(points)
    end
    if algorithm == "monotone_chain"
        return monotone_chain!(points)
    elseif algorithm == "monotone_chain_sorted"
        return monotone_chain!(points; sort=false)
    else
        error("the convex hull algorithm $algorithm is unknown")
    end
end

_zero_abs(x::Number) = x

_zero_abs(x::AbstractFloat) = x == -0.0 ? zero(x) : x

"""
    monotone_chain!(points::Vector{VN}; sort::Bool=true
                   ) where {N, VN<:AbstractVector{N}}

Compute the convex hull of a list of points in the plane using Andrew's
monotone-chain method.

### Input

- `points` -- list of 2D vectors; will be sorted in-place inside this method
- `sort`   -- (optional, default: `true`) flag for sorting the vertices
              lexicographically; sortedness is required for correctness

### Output

List of vectors containing the 2D coordinates of the corner points of the
convex hull.

### Notes

For large sets of points, it is convenient to use static vectors to get maximum
performance. For information on how to convert usual vectors into static
vectors, see the type `SVector` provided by the
[StaticArrays](http://juliaarrays.github.io/StaticArrays.jl/stable/)
package.

### Algorithm

This method implements Andrew's monotone-chain convex hull algorithm to
construct the convex hull of a set of ``n`` points in the plane in
``O(n \\log n)`` time. For further details see
[Monotone chain](https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain)
"""
function monotone_chain!(points::Vector{VN}; sort::Bool=true) where {N,VN<:AbstractVector{N}}
    @inline function build_hull!(semihull, iterator, points)
        @inbounds for i in iterator
            while length(semihull) >= 2 &&
                (right_turn(semihull[end - 1], semihull[end], points[i])
                 <=
                 zero(N))
                pop!(semihull)
            end
            push!(semihull, points[i])
        end
    end

    if sort
        # _zero_abs is necessary because floating-point arithmetic distinguishes
        # between 0.0 and -0.0, which leads to incorrect sorting (see #3455)
        # sort the points lexicographically
        sort!(points; by=x -> (_zero_abs(x[1]), _zero_abs(x[2])))
    end

    # build lower hull
    lower = Vector{VN}()
    iterator = eachindex(points)
    build_hull!(lower, iterator, points)

    # build upper hull
    upper = Vector{VN}()
    iterator = length(points):-1:1
    build_hull!(upper, iterator, points)

    # remove the last point of each segment because they are repeated
    copyto!(points, @view(lower[1:(end - 1)]))
    copyto!(points, length(lower), @view(upper[1:(end - 1)]))
    m = length(lower) + length(upper) - 2
    if m == 2 && _isapprox(points[1], points[2])
        # upper and lower chain consist of a single, identical point
        m = 1
    end
    return resize!(points, m)
end

@commutative function convex_hull(X::LazySet, ∅::EmptySet)
    return _convex_hull_emptyset(∅, X)
end

@commutative function convex_hull(X::LazySet, U::Universe)
    return _convex_hull_universe(U, X)
end

# ============== #
# disambiguation #
# ============== #

@commutative function convex_hull(U::Universe, ∅::EmptySet)
    return _convex_hull_universe(U, ∅)
end
