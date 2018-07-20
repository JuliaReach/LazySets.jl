import Base: ∈

export VPolygon

"""
    VPolygon{N<:Real} <: AbstractPolygon{N}

Type that represents a polygon by its vertices.

### Fields

- `vertices` -- the list of vertices

### Notes

The constructor of `VPolygon` runs a convex hull algorithm, and the given
vertices are sorted in counter-clockwise fashion.
The constructor flag `apply_convex_hull` can be used to skip the computation of
the convex hull.

- `VPolygon(vertices::Vector{Vector{N}};
            apply_convex_hull::Bool=true,
            algorithm::String="monotone_chain")`
"""
struct VPolygon{N<:Real} <: AbstractPolygon{N}
    vertices::Vector{Vector{N}}

    # default constructor that applies a convex hull algorithm
    function VPolygon(vertices::Vector{Vector{N}};
                      apply_convex_hull::Bool=true,
                      algorithm::String="monotone_chain") where {N<:Real}
        if apply_convex_hull
            return new{N}(convex_hull(vertices, algorithm=algorithm))
        else
            return new{N}(vertices)
        end
    end
end


# --- AbstractPolygon interface functions ---


"""
    tovrep(P::VPolygon{N})::VPolygon{N} where {N<:Real}

Build a vertex representation of the given polygon.

### Input

- `P` -- polygon in vertex representation

### Output

The identity, i.e., the same polygon instance.
"""
function tovrep(P::VPolygon{N})::VPolygon{N} where {N<:Real}
    return P
end

"""
    tohrep(P::VPolygon{N}, ::Type{HPOLYGON}=HPolygon
          )::AbstractHPolygon{N} where {N<:Real, HPOLYGON<:AbstractHPolygon}

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
               )::AbstractHPolygon{N} where {N<:Real, HPOLYGON<:AbstractHPolygon}
    n = length(vertices_list(P))
    if n == 0
        # no vertex -> no constraint
        constraints_list = Vector{LinearConstraint{N}}(0)
    elseif n == 1
        # only one vertex -> use function for singletons
        return convert(HPOLYGON, Singleton(P.vertices_list[1]))
    elseif n == 2
        # only two vertices -> use function for line segments
        return convert(HPOLYGON, LineSegment(P.vertices_list[1], P.vertices_list[2]))
    else
        # find right-most vertex
        i = div(n, 2)
        x = P.vertices_list[i][1]
        while i > 1 && P.vertices_list[i-1][1] > x
            # search forward in list
            i = i - 1
            x = P.vertices_list[i][1]
        end
        while i < n && P.vertices_list[i+1][1] > x
            # search backward in list
            i = i + 1
            x = P.vertices_list[i][1]
        end

        # create constraints counter-clockwise
        constraints_list = Vector{LinearConstraint{N}}(n)
        j = 1
        i_start = i
        i_stop = n
        for it in 1:2
            while i < i_stop
                constraints_list[j] =
                    halfspace_left(P.vertices_list[i], P.vertices_list[i+1])
                j += 1
                i += 1
            end
            if it == 1
                constraints_list[j] =
                    halfspace_left(P.vertices_list[n], P.vertices_list[1])
                j += 1
                i = 1
                i_stop = i_start
            end
        end
    end
    return HPOLYGON{N}(constraints_list)
end


# --- AbstractPolytope interface functions ---


"""
    vertices_list(P::VPolygon{N})::Vector{Vector{N}} where {N<:Real}

Return the list of vertices of a convex polygon in vertex representation.

### Input

- `P` -- a polygon vertex representation

### Output

List of vertices.
"""
function vertices_list(P::VPolygon{N})::Vector{Vector{N}} where {N<:Real}
    return P.vertices
end


# --- LazySet interface functions ---


"""
    σ(d::AbstractVector{<:Real}, P::VPolygon{N})::Vector{N} where {N<:Real}

Return the support vector of a polygon in a given direction.

### Input

- `d` -- direction
- `P` -- polygon in vertex representation

### Output

The support vector in the given direction.
If the direction has norm zero, the first vertex is returned.

### Algorithm

This implementation performs a brute-force search, comparing the projection of
each vector along the given direction.
It runs in ``O(n)`` where ``n`` is the number of vertices.

### Notes

For arbitrary points without structure this is the best one can do.
However, a more efficient approach can be used if the vertices of the polygon
have been sorted in counter-clockwise fashion.
In that case a binary search algorithm can be used that runs in ``O(\\log n)``.
See issue [#40](https://github.com/JuliaReach/LazySets.jl/issues/40).
"""
function σ(d::AbstractVector{<:Real},
           P::VPolygon{N})::Vector{N} where {N<:Real}
    if isempty(P.vertices)
        error("this polygon is empty")
    end
    i_max = 1
    @inbounds for i in 2:length(P.vertices)
        if dot(d, P.vertices[i] - P.vertices[i_max]) > zero(N)
            i_max = i
        end
    end
    return P.vertices[i_max]
end

"""
    an_element(P::VPolygon{N})::Vector{N} where {N<:Real}

Return some element of a polygon in vertex representation.

### Input

- `P` -- polygon in vertex representation

### Output

The first vertex of the polygon in vertex representation.
"""
function an_element(P::VPolygon{N})::Vector{N} where {N<:Real}
    if isempty(P.vertices)
        error("this polygon is empty")
    end
    return P.vertices[1]
end

"""
    ∈(x::AbstractVector{N}, P::VPolygon{N})::Bool where {N<:Real}

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

julia> ∈([4.5, 3.1], P)
false
julia> ∈([4.5, 3.0], P)
true
julia> ∈([4.4, 3.4], P)  #  point lies on the edge -> floating point error
false
julia> P = VPolygon([[2//1, 3//1], [3//1, 1//1], [5//1, 1//1], [4//1, 5//1]];
                     apply_convex_hull=false);

julia> ∈([44//10, 34//10], P)  #  with rational numbers the answer is correct
true
```
"""
function ∈(x::AbstractVector{N}, P::VPolygon{N})::Bool where {N<:Real}
    @assert length(x) == 2

    # special cases: 0 or 1 vertex
    if length(P.vertices) == 0
        return false
    elseif length(P.vertices) == 1
        return x == P.vertices[1]
    end

    zero_N = zero(N)
    if right_turn(P.vertices[1], x, P.vertices[end]) < zero_N
        return false
    end
    for i in 2:length(P.vertices)
        if right_turn(P.vertices[i], x, P.vertices[i-1]) < zero_N
            return false
        end
    end
    return true
end
