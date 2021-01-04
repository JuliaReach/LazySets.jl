import Base: rand,
             ∈

export VPolygon,
       remove_redundant_vertices,
       remove_redundant_vertices!,
       linear_map

"""
    VPolygon{N, VN<:AbstractVector{N}} <: AbstractPolygon{N}

Type that represents a polygon by its vertices.

### Fields

- `vertices` -- the list of vertices

### Notes

The constructor of `VPolygon` runs a convex hull algorithm on its vertices by
default, to remove the possibly redundant vertices. The vertices are sorted in
counter-clockwise fashion. Use the flag `apply_convex_hull=false` to skip the
computation of the convex hull.

- `VPolygon(vertices::Vector{Vector{N}};
            apply_convex_hull::Bool=true,
            algorithm::String="monotone_chain")`

### Examples

A polygon in vertex representation can be constructed by passing the list of
vertices. For example, we can build the right triangle

```jldoctest polygon_vrep
julia> P = VPolygon([[0, 0], [1, 0], [0, 1]]);

julia> P.vertices
3-element Array{Array{Int64,1},1}:
 [0, 0]
 [1, 0]
 [0, 1]
```

Alternatively, a `VPolygon` can be constructed passing a matrix of vertices,
where each *column* represents a vertex:

```jldoctest polygon_vrep
julia> M = [0 1 0; 0 0 1.]
2×3 Array{Float64,2}:
 0.0  1.0  0.0
 0.0  0.0  1.0

julia> P = VPolygon(M);

julia> P.vertices
3-element Array{Array{Float64,1},1}:
 [0.0, 0.0]
 [1.0, 0.0]
 [0.0, 1.0]
```
"""
struct VPolygon{N, VN<:AbstractVector{N}} <: AbstractPolygon{N}
    vertices::Vector{VN}

    # default constructor that applies a convex hull algorithm
    function VPolygon(vertices::Vector{VN};
                      apply_convex_hull::Bool=true,
                      algorithm::String="monotone_chain"
                     ) where {N, VN<:AbstractVector{N}}
        if apply_convex_hull
            return new{N, VN}(convex_hull(vertices, algorithm=algorithm))
        else
            return new{N, VN}(vertices)
        end
    end
end

isoperationtype(::Type{<:VPolygon}) = false
isconvextype(::Type{<:VPolygon}) = true

# constructor with empty vertices list
VPolygon{N}() where {N} = VPolygon(Vector{Vector{N}}(), apply_convex_hull=false)

# constructor with no vertices of type Float64
VPolygon() = VPolygon{Float64}()

# constructor from rectangular matrix
function VPolygon(vertices_matrix::MT; apply_convex_hull::Bool=true,
                  algorithm::String="monotone_chain") where {MT<:AbstractMatrix}
    @assert size(vertices_matrix, 1) == 2 "the number of rows of the matrix of vertices " *
                                           "should be 2, but it is $(size(vertices_matrix, 1))"

    vertices = [vertices_matrix[:, j] for j in 1:size(vertices_matrix, 2)]
    return VPolygon(vertices; apply_convex_hull=apply_convex_hull, algorithm=algorithm)
end

"""
    remove_redundant_vertices!(P::VPolygon;
                               [algorithm]::String="monotone_chain")

Remove the redundant vertices of the given polygon.

### Input

- `P`         -- polygon in vertex representation
- `algorithm` -- (optional, default: "monotone_chain") the algorithm used to
                 compute the convex hull

### Output

A new polygon such that its vertices are the convex hull of the given polygon.

### Algorithm

A convex hull algorithm is used to compute the convex hull of the vertices of the
given input polygon `P`; see `?convex_hull` for details on the available algorithms.
The vertices of the output polygon are sorted in counter-clockwise fashion.
"""
function remove_redundant_vertices!(P::VPolygon;
                                    algorithm::String="monotone_chain")
    convex_hull!(P.vertices; algorithm=algorithm)
end

"""
    remove_redundant_vertices(P::VPolygon;
                              [algorithm]::String="monotone_chain")

Return the polygon obtained by removing the redundant vertices of the given polygon.

### Input

- `P`         -- polygon in vertex representation
- `algorithm` -- (optional, default: "monotone_chain") the algorithm used to
                 compute the convex hull

### Output

A new polygon such that its vertices are the convex hull of the given polygon.

### Algorithm

A convex hull algorithm is used to compute the convex hull of the vertices of the
given input polygon `P`; see `?convex_hull` for details on the available algorithms.
The vertices of the output polygon are sorted in counter-clockwise fashion.
"""
function remove_redundant_vertices(P::VPolygon;
                                   algorithm::String="monotone_chain")
    return remove_redundant_vertices!(copy(P), algorithm=algorithm)
end

function linear_map(M::AbstractMatrix, P::VPolygon)
    @assert size(M, 2) == 2 "a linear map of size $(size(M)) cannot be applied to a set of dimension 2"
    return _linear_map_vrep(M, P)
end

@inline function _linear_map_vrep(M::AbstractMatrix, P::VPolygon)
    return broadcast(v -> M * v, vertices_list(P)) |> VPolygon
end

# --- AbstractPolygon interface functions ---


"""
    tovrep(P::VPolygon)

Build a vertex representation of the given polygon.

### Input

- `P` -- polygon in vertex representation

### Output

The identity, i.e., the same polygon instance.
"""
function tovrep(P::VPolygon)
    return P
end

"""
    tohrep(P::VPolygon{N}, ::Type{HPOLYGON}=HPolygon
          ) where {N, HPOLYGON<:AbstractHPolygon}

Build a constraint representation of the given polygon.

### Input

- `P`        -- polygon in vertex representation
- `HPOLYGON` -- (optional, default: `HPolygon`) type of target polygon

### Output

The same polygon but in constraint representation, an `AbstractHPolygon`.

### Algorithm

The algorithms consists of adding an edge for each consecutive pair of vertices.
Since the vertices are already ordered in counter-clockwise fashion (CWW), the
constraints will be sorted automatically (CCW) if we start with the first edge
between the first and second vertex.
"""
function tohrep(P::VPolygon{N}, ::Type{HPOLYGON}=HPolygon
               ) where {N, HPOLYGON<:AbstractHPolygon}
    vl = vertices_list(P)
    n = length(vl)
    if n == 0
        # no vertex -> empy set
        return EmptySet{N}(2)
    elseif n == 1
        # only one vertex -> use function for singletons
        return convert(HPOLYGON, Singleton(vl[1]))
    elseif n == 2
        # only two vertices -> use function for line segments
        return convert(HPOLYGON, LineSegment(vl[1], vl[2]))
    else
        # find right-most vertex
        i = div(n, 2)
        x = vl[i][1]
        while i > 1 && vl[i-1][1] > x
            # search forward in list
            i = i - 1
            x = vl[i][1]
        end
        while i < n && vl[i+1][1] > x
            # search backward in list
            i = i + 1
            x = vl[i][1]
        end

        # create constraints ordered in CCW starting at the right-most index
        upper_hull = [halfspace_left(vl[j], vl[j+1]) for j in i:length(vl)-1]
        mid_hull = [halfspace_left(vl[end], vl[1])]
        lower_hull = [halfspace_left(vl[j], vl[j+1]) for j in 1:i-1]
        constraints_list = vcat(upper_hull, mid_hull, lower_hull)
    end
    return HPOLYGON(constraints_list)
end


# --- AbstractPolytope interface functions ---


"""
    vertices_list(P::VPolygon)

Return the list of vertices of a convex polygon in vertex representation.

### Input

- `P` -- a polygon vertex representation

### Output

List of vertices.
"""
function vertices_list(P::VPolygon)
    return P.vertices
end


# --- LazySet interface functions ---


"""
    σ(d::AbstractVector, P::VPolygon)

Return the support vector of a polygon in a given direction.

### Input

- `d` -- direction
- `P` -- polygon in vertex representation

### Output

The support vector in the given direction.
If the direction has norm zero, the first vertex is returned.

### Algorithm

This implementation uses a binary search algorithm when the polygon has more
than ten vertices and a brute-force search when it has ten or less.
For the brute-force search, it compares the projection of
each vector along the given direction and runs in ``O(n)`` where
``n`` is the number of vertices.
For the binary search the algorithm runs in ``O(log n)``.
We follow [this implementation](http://geomalgorithms.com/a14-_extreme_pts.html#polyMax_2D())
based on an algorithm described in [1].

[1] Joseph O'Rourke, Computational Geometry in C (2nd Edition)
"""
function σ(d::AbstractVector, P::VPolygon)
    @assert !isempty(P.vertices) "the polygon has no vertices"
    return P.vertices[_σ_helper(d, P)]
end

# heuristic to define the method used to compute the support vector of a polygon
# in vertex representation; if the number of vertices of the polygon is smaller
# than this value, the brute force method is used; otherwise binary search is used
const binary_or_brute_force = 10

# return the index of a support vector of P along d from the list of vertices of P
function _σ_helper(d::AbstractVector, P::VPolygon)
    vlist = P.vertices
    return _σ_helper(d, vlist)
end

function _σ_helper(d::AbstractVector, vlist::Vector{VN}) where {N, VN<:AbstractVector{N}}
    if length(vlist) > binary_or_brute_force
        return _binary_support_vector(d, vlist)
    else
        return _brute_force_support_vector(d, vlist)
    end
end

function _brute_force_support_vector(d::AbstractVector, P::VPolygon)
    return _brute_force_support_vector(d, P.vertices)
end

function _brute_force_support_vector(d::AbstractVector{M}, vlist::Vector{VT}) where {M, T, VT<:AbstractVector{T}}
    i_max = 1
    N = promote_type(M, T)
    @inbounds for i in 2:length(vlist)
        if dot(d, vlist[i] - vlist[i_max]) > zero(N)
            i_max = i
        end
    end
    return i_max
end

function _binary_support_vector(d::AbstractVector, P::VPolygon)
    return _binary_support_vector(d, P.vertices)
end

function _binary_support_vector(d::AbstractVector, vlist::Vector{VT}) where {T, VT<:AbstractVector{T}}
    m = length(vlist)
    @assert m > 2 "the number of vertices in the binary support vector approach should " *
                  "be at least three, but it is $m"

    # add extra vertex on the end equal to the first
    push!(vlist, vlist[1])

    # start chain = [1,n+1] with P.vertices[n+1] = P.vertices[1]
    a = 1; b = m + 1
    A = vlist[2] - vlist[1]
    upA = _up(d, A)

    # test if P.vertices[0] is a local maximum
    if (!upA && !_above(d, vlist[m], vlist[1]))
        # if vlist[1] is the maximum, remove the extra point added
        pop!(vlist)
        return 1
    end
    while true
        # midpoint of [a,b], and 1<c<n+1
        c = round(Int, (a + b) / 2)
        C = vlist[c + 1] - vlist[c]
        upC = _up(d, C)
        if (!upC && !_above(d, vlist[c - 1], vlist[c]))
            # vlist[c] is a local maximum, remove the extra point added
            pop!(vlist)
            return c
        end

        # no max yet, so continue with the binary search
        # pick one of the two subchains [a,c] or [c,b]
        if (upA && upC && !_above(d, vlist[a], vlist[c])) ||
        (!upA && (upC || (!upC && _above(d, vlist[a], vlist[c]))))
            a = c
            A = C
            upA = upC
        else
            b = c
        end
        if (b <= a + 1) # the chain is impossibly small
            pop!(vlist) # remove the extra point added
            throw(ErrorException("something went wrong")) # return an error
        end
    end
end

"""
    an_element(P::VPolygon)

Return some element of a polygon in vertex representation.

### Input

- `P` -- polygon in vertex representation

### Output

The first vertex of the polygon in vertex representation.
"""
function an_element(P::VPolygon)
    @assert !isempty(P.vertices) "the polygon has no vertices"
    return P.vertices[1]
end

"""
    ∈(x::AbstractVector, P::VPolygon)

Check whether a given point is contained in a polygon in vertex representation.

### Input

- `x` -- point/vector
- `P` -- polygon in vertex representation

### Output

`true` iff ``x ∈ P``.

### Algorithm

This implementation exploits that the polygon's vertices are sorted in
counter-clockwise fashion.
Under this assumption we can just check if the vertex lies on the left of each
edge, using the dot product.

### Examples

```jldoctest
julia> P = VPolygon([[2.0, 3.0], [3.0, 1.0], [5.0, 1.0], [4.0, 5.0]];
                    apply_convex_hull=false);

julia> [4.5, 3.1] ∈ P
false
julia> [4.5, 3.0] ∈ P
true
julia> [4.4, 3.4] ∈ P  #  point lies on the edge
true
```
"""
function ∈(x::AbstractVector, P::VPolygon)
    @assert length(x) == 2

    # special cases: 0 or 1 vertex
    if length(P.vertices) == 0
        return false
    elseif length(P.vertices) == 1
        return x == P.vertices[1]
    end

    if !is_right_turn(P.vertices[1], x, P.vertices[end])
        return false
    end
    for i in 2:length(P.vertices)
        if !is_right_turn(P.vertices[i], x, P.vertices[i-1])
            return false
        end
    end
    return true
end

"""
    _random_zero_sum_vector(rng::AbstractRNG, N::Type{<:Real}, n::Int)

Create a random vector with entries whose sum is zero.

### Input

- `rng` -- random number generator
- `N`   -- numeric type
- `n`   -- length of vector

### Output

A random vector of random numbers such that all positive entries come first and
all negative entries come last, and such that the total sum is zero.

### Algorithm

This is a preprocessing step of the algorithm
[here](https://stackoverflow.com/a/47358689) based on
[P. Valtr. Probability that n random points are in convex
position](https://link.springer.com/article/10.1007%2FBF01271274).
"""
function _random_zero_sum_vector(rng::AbstractRNG, N::Type{<:Real}, n::Int)
    # generate a sorted list of random x and y coordinates
    list = sort!(randn(rng, N, n))
    while (length(remove_duplicates_sorted!(list)) < n)
        # make sure that no duplicates exist
        list = sort!(append!(list, randn(rng, N, length(list) - n)))
    end
    # lists of consecutive points
    in_l1 = rand(rng, Bool, n-2)
    l1 = Vector{N}() # normal
    l2 = Vector{N}() # inverted
    push!(l1, list[1])
    push!(l2, list[1])
    for i in 2:n-1
        push!(in_l1[i-1] ? l1 : l2, list[i])
    end
    push!(l1, list[end])
    push!(l2, list[end])
    # convert to vectors representing the distance (order does not matter)
    dist = Vector{N}()
    sizehint!(dist, n)
    for i in 1:length(l1)-1
        push!(dist, l1[i+1] - l1[i])
    end
    for i in 1:length(l2)-1
        push!(dist, l2[i] - l2[i+1])
    end
    @assert isapprox(sum(dist), zero(N), atol=1e-6)
    return dist
end

"""
    rand(::Type{VPolygon}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

Create a random polygon in vertex representation.

### Input

- `VPolygon`     -- type for dispatch
- `N`            -- (optional, default: `Float64`) numeric type
- `dim`          -- (optional, default: 2) dimension
- `rng`          -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`         -- (optional, default: `nothing`) seed for reseeding
- `num_vertices` -- (optional, default: `-1`) number of vertices of the
                    polygon (see comment below)

### Output

A random polygon in vertex representation.

### Algorithm

We follow the idea [here](https://stackoverflow.com/a/47358689) based on
[P. Valtr. Probability that n random points are in convex
position](https://link.springer.com/article/10.1007%2FBF01271274).
There is also a nice video available
[here](http://cglab.ca/~sander/misc/ConvexGeneration/convex.html).

The number of vertices can be controlled with the argument `num_vertices`.
For a negative value we choose a random number in the range `3:10`.
"""
function rand(::Type{VPolygon};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int, Nothing}=nothing,
              num_vertices::Int=-1)
    @assert dim == 2 "cannot create a random VPolygon of dimension $dim"
    rng = reseed(rng, seed)
    if num_vertices < 0
        num_vertices = rand(3:10)
    end

    # special cases, 0 or 1 vertex
    if num_vertices == 0
        return VPolygon{N}()
    elseif num_vertices == 1
        return VPolygon([randn(rng, N, 2)])
    end

    # general case, >= 2 vertices

    # get random horizontal and vertical vectors
    horiz = _random_zero_sum_vector(rng, N, num_vertices)
    vert = _random_zero_sum_vector(rng, N, num_vertices)

    # randomly combine horizontal and vertical vectors
    m = num_vertices
    directions = Vector{Vector{N}}(undef, num_vertices)
    shuffle(rng, vert)
    for (i, x) in enumerate(horiz)
        y = splice!(vert, rand(rng, 1:m))
        directions[i] = [x, y]
        m -= 1
    end
    sort!(directions, lt=<=) # sort by angle

    # connect directions
    vertices = Vector{Vector{N}}(undef, num_vertices)
    # random starting point
    vertices[1] = randn(rng, N, 2)
    for i in 1:length(directions)-1
        vertices[i+1] = vertices[i] + directions[i]
    end
    @assert isapprox(vertices[end] + directions[end], vertices[1], atol=1e-6)
    return VPolygon(vertices; apply_convex_hull=true)
end

"""
    constraints_list(P::VPolygon)

Return the list of constraints defining a polygon in V-representation.

### Input

- `P` -- polygon in V-representation

### Output

The list of constraints of the polygon.

### Algorithm

First the H-representation of ``P`` is computed, then its list of constraints
is returned.
"""
function constraints_list(P::VPolygon)
    return constraints_list(tohrep(P))
end

"""
    translate(P::VPolygon, v::AbstractVector)

Translate (i.e., shift) a polygon in vertex representation by a given vector.

### Input

- `P` -- polygon in vertex representation
- `v` -- translation vector

### Output

A translated polygon in vertex representation.

### Algorithm

We add the vector to each vertex of the polygon.
"""
function translate(P::VPolygon, v::AbstractVector)
    @assert length(v) == dim(P) "cannot translate a $(dim(P))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    return VPolygon([x + v for x in vertices_list(P)])
end

# =======================
# Concrete intersection
# =======================

function _isinside(p, a, b)
    @inbounds begin
        α = (b[1] - a[1]) * (p[2] - a[2])
        β = (b[2] - a[2]) * (p[1] - a[1])
    end
    return !_leq(α, β)
end

function _intersection_line_segments(a, b, s, f)
    @inbounds begin
        dc = [a[1] - b[1], a[2] - b[2]]
        dp = [s[1] - f[1], s[2] - f[2]]
        n1 = a[1] * b[2] - a[2] * b[1]
        n2 = s[1] * f[2] - s[2] * f[1]
        n3 = 1.0 / (dc[1] * dp[2] - dc[2] * dp[1])

        α = (n1 * dp[1] - n2 * dc[1]) * n3
        β = (n1 * dp[2] - n2 * dc[2]) * n3
    end
    return [α, β]
end

function _intersection_vrep(spoly::AbstractVector{VT},
                            cpoly::AbstractVector{VT}) where {VT<:AbstractVector}
    outarr = spoly
    q = cpoly[end]
    for p in cpoly
        inarr = outarr
        outarr = Vector{VT}()
        isempty(inarr) && break
        s = inarr[end]
        for vtx in inarr
            if _isinside(vtx, q, p)
                if !_isinside(s, q, p)
                    push!(outarr, _intersection_line_segments(q, p, s, vtx))
                end
                push!(outarr, vtx)
            elseif _isinside(s, q, p)
                push!(outarr, _intersection_line_segments(q, p, s, vtx))
            end
            s = vtx
        end
        q = p
    end
    return outarr
end
