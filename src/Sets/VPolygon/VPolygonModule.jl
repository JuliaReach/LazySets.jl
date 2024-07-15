module VPolygonModule

using Reexport, Requires

using ..LazySets: AbstractPolygon, AbstractHPolygon, HPolygon, halfspace_left,
                  is_right_turn, _area_vlist, _linear_map_vrep
using Random: AbstractRNG, GLOBAL_RNG, shuffle
using ReachabilityBase.Arrays: isabove, rand_pos_neg_zerosum_vector
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Require: require
using LinearAlgebra: dot

@reexport import ..API: an_element, area, constraints_list, isoperationtype,
                        rand, vertices_list, ∈, linear_map, permute, project,
                        σ, translate, translate!
@reexport import ..LazySets: remove_redundant_vertices,
                             remove_redundant_vertices!, tohrep, tovrep
@reexport using ..API

export VPolygon

# heuristic to define the method used to compute the support vector of a polygon
# in vertex representation; if the number of vertices of the polygon is smaller
# than this value, the brute force method is used; otherwise binary search is used
const BINARY_OR_BRUTE_FORCE = 10

# range of the default number of vertices in `rand`
const DEFAULT_RAND_VERTEX_RANGE = 3:10

"""
    VPolygon{N, VN<:AbstractVector{N}} <: AbstractPolygon{N}

Type that represents a polygon by its vertices.

### Fields

- `vertices` -- the list of vertices

### Notes

This type assumes that all vertices are sorted in counter-clockwise fashion.

To ensure this property, the constructor of `VPolygon` runs a convex-hull
algorithm on the vertices by default. This also removes redundant vertices.
If the vertices are known to be sorted, the flag `apply_convex_hull=false` can
be used to skip this preprocessing.

### Examples

A polygon in vertex representation can be constructed by passing the list of
vertices. For example, we can build the right triangle

```jldoctest polygon_vrep
julia> P = VPolygon([[0, 0], [1, 0], [0, 1]]);

julia> P.vertices
3-element Vector{Vector{Int64}}:
 [0, 0]
 [1, 0]
 [0, 1]
```

Alternatively, a `VPolygon` can be constructed passing a matrix of vertices,
where each *column* represents a vertex:

```jldoctest polygon_vrep
julia> M = [0 1 0; 0 0 1.]
2×3 Matrix{Float64}:
 0.0  1.0  0.0
 0.0  0.0  1.0

julia> P = VPolygon(M);

julia> P.vertices
3-element Vector{Vector{Float64}}:
 [0.0, 0.0]
 [1.0, 0.0]
 [0.0, 1.0]
```
"""
struct VPolygon{N,VN<:AbstractVector{N}} <: AbstractPolygon{N}
    vertices::Vector{VN}

    # default constructor that applies a convex hull algorithm
    function VPolygon(vertices::Vector{VN};
                      apply_convex_hull::Bool=true,
                      algorithm::String="monotone_chain") where {N,VN<:AbstractVector{N}}
        if apply_convex_hull
            vertices = convex_hull(vertices; algorithm=algorithm)
        end
        return new{N,VN}(vertices)
    end
end

isoperationtype(::Type{<:VPolygon}) = false

# constructor with empty vertices list
VPolygon{N}() where {N} = VPolygon(Vector{Vector{N}}(); apply_convex_hull=false)

# constructor with no vertices of type Float64
VPolygon() = VPolygon{Float64}()

# constructor from rectangular matrix
function VPolygon(vertices_matrix::MT; apply_convex_hull::Bool=true,
                  algorithm::String="monotone_chain") where {MT<:AbstractMatrix}
    @assert size(vertices_matrix, 1) == 2 "the number of rows of the matrix " *
                                          "of vertices should be 2, but it is $(size(vertices_matrix, 1))"

    vertices = [vertices_matrix[:, j] for j in axes(vertices_matrix, 2)]
    return VPolygon(vertices; apply_convex_hull=apply_convex_hull,
                    algorithm=algorithm)
end

"""
    remove_redundant_vertices!(P::VPolygon;
                               [algorithm]::String="monotone_chain")

Remove the redundant vertices from the given polygon in-place.

### Input

- `P`         -- polygon in vertex representation
- `algorithm` -- (optional, default: "monotone_chain") the algorithm used to
                 compute the convex hull

### Output

The modified polygon whose redundant vertices have been removed.

### Algorithm

A convex-hull algorithm is used to compute the convex hull of the vertices of
the polygon `P`; see `?convex_hull` for details on the available algorithms.
The vertices are sorted in counter-clockwise fashion.
"""
function remove_redundant_vertices!(P::VPolygon;
                                    algorithm::String="monotone_chain")
    convex_hull!(P.vertices; algorithm=algorithm)
    return P
end

"""
    remove_redundant_vertices(P::VPolygon;
                              [algorithm]::String="monotone_chain")

Return a polygon obtained by removing the redundant vertices of the given
polygon.

### Input

- `P`         -- polygon in vertex representation
- `algorithm` -- (optional, default: "monotone_chain") the algorithm used to
                 compute the convex hull

### Output

A new polygon such that its vertices are the convex hull of the given polygon.

### Algorithm

See [`remove_redundant_vertices!(::VPolygon)`](@ref).
"""
function remove_redundant_vertices(P::VPolygon;
                                   algorithm::String="monotone_chain")
    return remove_redundant_vertices!(copy(P); algorithm=algorithm)
end

"""
    linear_map(M::AbstractMatrix, P::VPolygon; [apply_convex_hull]::Bool=false)

Concrete linear map of a polygon in vertex representation.

### Input

- `M`                 -- matrix
- `P`                 -- polygon in vertex representation
- `apply_convex_hull` -- (optional; default: `false`) flag to apply a
                         convex-hull operation (only relevant for
                         higher-dimensional maps)

### Output

The type of the result depends on the dimension. in 1D it is an interval, in 2D
it is a `VPolygon`, and in all other cases it is a `VPolytope`.

### Algorithm

This implementation uses the internal `_linear_map_vrep` method.
"""
function linear_map(M::AbstractMatrix, P::VPolygon;
                    apply_convex_hull::Bool=false)
    @assert size(M, 2) == 2 "a linear map of size $(size(M)) cannot be " *
                            "applied to a set of dimension 2"
    return _linear_map_vrep(M, P; apply_convex_hull=apply_convex_hull)
end

"""
    tovrep(P::VPolygon)

Build a vertex representation of the given polygon.

### Input

- `P` -- polygon in vertex representation

### Output

The same polygon instance.
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

A polygon in constraint representation, an `AbstractHPolygon`.

### Algorithm

The algorithm adds an edge for each consecutive pair of vertices.
Since the vertices are already ordered in counter-clockwise fashion (CCW), the
constraints will be sorted automatically (CCW).
"""
function tohrep(P::VPolygon{N}, ::Type{HPOLYGON}=HPolygon) where {N,HPOLYGON<:AbstractHPolygon}
    vl = P.vertices
    n = length(vl)
    if n == 0
        # no vertex
        return EmptySet{N}(2)
    elseif n == 1
        # only one vertex -> use function for singletons
        require(@__MODULE__, :LauySets; fun_name="convert")
        return convert(HPOLYGON, Singleton(vl[1]))
    elseif n == 2
        # only two vertices -> use function for line segments
        require(@__MODULE__, :LauySets; fun_name="convert")
        return convert(HPOLYGON, LineSegment(vl[1], vl[2]))
    else
        # find right-most vertex
        i = div(n, 2)
        x = vl[i][1]
        while i > 1 && vl[i - 1][1] > x
            # search forward in list
            i = i - 1
            x = vl[i][1]
        end
        while i < n && vl[i + 1][1] > x
            # search backward in list
            i = i + 1
            x = vl[i][1]
        end

        # create constraints ordered in CCW starting at the right-most index
        upper_hull = [halfspace_left(vl[j], vl[j + 1]) for j in i:(length(vl) - 1)]
        mid_hull = [halfspace_left(vl[end], vl[1])]
        lower_hull = [halfspace_left(vl[j], vl[j + 1]) for j in 1:(i - 1)]
        constraints_list = vcat(upper_hull, mid_hull, lower_hull)
    end
    return HPOLYGON(constraints_list)
end

"""
    vertices_list(P::VPolygon; kwargs...)

Return the list of vertices of a polygon in vertex representation.

### Input

- `P` -- polygon in vertex representation

### Output

The list of vertices.
"""
function vertices_list(P::VPolygon; kwargs...)
    return P.vertices
end

"""
    σ(d::AbstractVector, P::VPolygon)

Return a support vector of a polygon in a given direction.

### Input

- `d` -- direction
- `P` -- polygon in vertex representation

### Output

A support vector in the given direction.
If the direction has norm zero, the first vertex is returned.

### Algorithm

This implementation uses a binary search algorithm when the polygon has more
than $BINARY_OR_BRUTE_FORCE vertices and a brute-force search when it has
$BINARY_OR_BRUTE_FORCE or fewer vertices.
The brute-force search compares the projection of each vector along the given
direction and runs in ``O(n)`` where ``n`` is the number of vertices.
The binary search runs in ``O(log n)`` and we follow
[this implementation](http://geomalgorithms.com/a14-_extreme_pts.html#polyMax_2D())
based on an algorithm described in [1].

[1] Joseph O'Rourke, *Computational Geometry in C (2nd Edition)*.
"""
function σ(d::AbstractVector, P::VPolygon)
    @assert !isempty(P.vertices) "the polygon has no vertices"
    return P.vertices[_σ_helper(d, P)]
end

# return the index of a support vector of P along d from the vertices list of P
function _σ_helper(d::AbstractVector, P::VPolygon)
    vlist = P.vertices
    return _σ_helper(d, vlist)
end

function _σ_helper(d::AbstractVector, vlist::Vector{VN}) where {N,VN<:AbstractVector{N}}
    if length(vlist) > BINARY_OR_BRUTE_FORCE
        return _binary_support_vector(d, vlist)
    else
        return _brute_force_support_vector(d, vlist)
    end
end

function _brute_force_support_vector(d::AbstractVector, P::VPolygon)
    return _brute_force_support_vector(d, P.vertices)
end

function _brute_force_support_vector(d::AbstractVector{M},
                                     vlist::Vector{VT}) where {M,T,VT<:AbstractVector{T}}
    max_idx = 1
    @inbounds max_ρ = dot(d, vlist[1])
    @inbounds for i in 2:length(vlist)
        ρ_i = dot(d, vlist[i])
        if ρ_i > max_ρ
            max_idx = i
            max_ρ = ρ_i
        end
    end
    return max_idx
end

function _binary_support_vector(d::AbstractVector, P::VPolygon)
    return _binary_support_vector(d, P.vertices)
end

# checks if the given vector is pointing toward the given direction
@inline function _similar_direction(u::AbstractVector, v::AbstractVector)
    return dot(u, v) > 0
end

function _binary_support_vector(d::AbstractVector,
                                vlist::Vector{VT}) where {T,VT<:AbstractVector{T}}
    m = length(vlist)
    @assert m > 2 "the number of vertices in the binary-search approach " *
                  "should be at least three, but it is $m"

    @inbounds begin
        # add the first vertex at the end again
        push!(vlist, vlist[1])

        # start chain = [1,n+1] with vlist[n+1] = vlist[1]
        a = 1
        b = m + 1
        A = vlist[2] - vlist[1]
        upA = _similar_direction(d, A)
        if (!upA && !isabove(d, vlist[m], vlist[1]))
            # vlist[1] is the maximum
            pop!(vlist)  # remove the extra point added
            return 1
        end
        while true
            # midpoint of [a,b], and 1<c<n+1
            c = round(Int, (a + b) / 2)
            C = vlist[c + 1] - vlist[c]
            upC = _similar_direction(d, C)
            if (!upC && !isabove(d, vlist[c - 1], vlist[c]))
                # vlist[c] is the maximum
                pop!(vlist)  # remove the extra point added
                return c
            end

            # no max yet, so continue with the binary search
            # pick one of the two subchains [a,c] or [c,b]
            if (upA && upC && !isabove(d, vlist[a], vlist[c])) ||
               (!upA && (upC || (!upC && isabove(d, vlist[a], vlist[c]))))
                a = c
                A = C
                upA = upC
            else
                b = c
            end
            if (b <= a + 1)  # the chain is impossibly small
                pop!(vlist)  # remove the extra point added
                error("something went wrong")
            end
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
julia> P = VPolygon([[2.0, 3.0], [3.0, 1.0], [5.0, 1.0], [4.0, 5.0]]);

julia> [4.5, 3.1] ∈ P
false
julia> [4.5, 3.0] ∈ P
true
julia> [4.4, 3.4] ∈ P  #  point lies on the edge
true
```
"""
function ∈(x::AbstractVector, P::VPolygon)
    @assert length(x) == 2 "a $(length(x))-dimensional vector is " *
                           "incompatible with an $(dim(P))-dimensional set"

    # special cases: 0 or 1 vertex
    @inbounds begin
        if length(P.vertices) == 0
            return false
        elseif length(P.vertices) == 1
            return x == P.vertices[1]
        end

        if !is_right_turn(P.vertices[1], x, P.vertices[end])
            return false
        end
        for i in 2:length(P.vertices)
            if !is_right_turn(P.vertices[i], x, P.vertices[i - 1])
                return false
            end
        end
    end
    return true
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

### Notes

The number of vertices can be controlled with the argument `num_vertices`.
For a negative value we choose a random number in the range
`$DEFAULT_RAND_VERTEX_RANGE`.

### Algorithm

We follow the idea described [here](https://stackoverflow.com/a/47358689) based
on [1]. There is also a nice video available
[here](http://cglab.ca/~sander/misc/ConvexGeneration/convex.html).

[1] Pavel Valtr: *Probability that n random points are in convex position*.
Discret. Comput. Geom. 1995.
"""
function rand(::Type{VPolygon};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing,
              num_vertices::Int=-1)
    @assert dim == 2 "cannot create a random VPolygon of dimension $dim"
    rng = reseed!(rng, seed)
    if num_vertices < 0
        num_vertices = rand(rng, DEFAULT_RAND_VERTEX_RANGE)
    end

    # special cases, 0 or 1 vertex
    if num_vertices == 0
        return VPolygon{N}()
    elseif num_vertices == 1
        return VPolygon([randn(rng, N, 2)])
    end

    # general case, >= 2 vertices

    # get random horizontal and vertical vectors
    horiz = rand_pos_neg_zerosum_vector(num_vertices; N=N, rng=rng)
    vert = rand_pos_neg_zerosum_vector(num_vertices; N=N, rng=rng)

    # randomly combine horizontal and vertical vectors
    m = num_vertices
    directions = Vector{Vector{N}}(undef, num_vertices)
    shuffle(rng, vert)
    for (i, x) in enumerate(horiz)
        y = splice!(vert, rand(rng, 1:m))
        directions[i] = [x, y]
        m -= 1
    end
    sort!(directions; lt=<=) # sort by angle

    # connect directions
    vertices = Vector{Vector{N}}(undef, num_vertices)
    # random starting point
    @inbounds begin
        vertices[1] = randn(rng, N, 2)
        for i in 1:(length(directions) - 1)
            vertices[i + 1] = vertices[i] + directions[i]
        end
        @assert isapprox(vertices[end] + directions[end], vertices[1], atol=1e-6)
    end
    return VPolygon(vertices; apply_convex_hull=true)
end

"""
    constraints_list(P::VPolygon)

Return the list of constraints defining a polygon in vertex representation.

### Input

- `P` -- polygon in vertex representation

### Output

The list of constraints of the polygon.

### Algorithm

We convert to constraint representation using `tohrep`.
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

### Notes

See [`translate!(::VPolygon, ::AbstractVector)`](@ref) for the in-place version.
"""
function translate(P::VPolygon, v::AbstractVector)
    return translate!(deepcopy(P), v)
end

"""
    translate!(P::VPolygon, v::AbstractVector)

Translate (i.e., shift) a polygon in vertex representation by a given vector,
in-place.

### Input

- `P` -- polygon in vertex representation
- `v` -- translation vector

### Output

The polygon translated by the vector.

### Algorithm

We add the vector to each vertex of the polygon.

### Notes

See [`translate(::VPolygon, ::AbstractVector)`](@ref) for the out-of-place
version.
"""
function translate!(P::VPolygon, v::AbstractVector)
    @assert length(v) == dim(P) "cannot translate a $(dim(P))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    for x in P.vertices
        x .+= v
    end
    return P
end

function project(V::VPolygon, block::AbstractVector{Int}; kwargs...)
    if length(block) == 1
        @assert block[1] == 1 || block[1] == 2 "invalid projection to $block"
        l, h = extrema(V, block[1])
        return Interval(l, h)
    elseif length(block) == 2
        if block[1] == 1 && block[2] == 2
            return V  # no projection
        else
            @assert block[1] == 2 && block[2] == 1 "invalid projection to $block"
            return permute(V, block)  # swap dimensions
        end
    end
    throw(ArgumentError("invalid projection to $block"))
end

"""
    permute(V::VPolygon, p::AbstractVector{Int})

Permute the dimensions according to a permutation vector.

### Input

- `P` -- polygon in vertex representation
- `p` -- permutation vector

### Output

The permuted polygon in vertex representation.
"""
function permute(V::VPolygon, p::AbstractVector{Int})
    return VPolygon([v[p] for v in V.vertices]; apply_convex_hull=true)
end

"""
    area(V::VPolygon)

Compute the area of a polygon in vertex representation.

### Input

- `V` -- polygon in vertex representation

### Output

A number representing the area of `V`.

### Algorithm

See [`area(::LazySet)`](@ref).
"""
function area(V::VPolygon)
    return _area_vlist(V.vertices; apply_convex_hull=false)
end

include("init.jl")

end  # module
