import Base: ∈

export VPolygon

"""
    VPolygon{N<:Real} <: AbstractPolygon{N}

Type that represents a polygon by its vertices.

### Fields

- `vertices_list` -- the list of vertices

### Notes

The constructor of `VPolygon` runs a convex hull algorithm, and the given
vertices are sorted in counter-clockwise fashion.
The constructor flag `apply_convex_hull` can be used to skip the computation of
the convex hull.

- `VPolygon(vertices_list::Vector{Vector{N}};
            apply_convex_hull::Bool=true,
            algorithm::String="monotone_chain")`
"""
struct VPolygon{N<:Real} <: AbstractPolygon{N}
    vertices_list::Vector{Vector{N}}

    # default constructor that applies a convex hull algorithm
    function VPolygon(vertices_list::Vector{Vector{N}};
                      apply_convex_hull::Bool=true,
                      algorithm::String="monotone_chain") where {N<:Real}
        if apply_convex_hull
            return new{N}(convex_hull(vertices_list, algorithm=algorithm))
        else
            return new{N}(vertices_list)
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
    tohrep(P::VPolygon{N})::AbstractHPolygon{N} where {N<:Real}

Build a constraint representation of the given polygon.

### Input

- `P` -- polygon in vertex representation

### Output

The same polygon but in constraint representation, an `AbstractHPolygon`.
"""
function tohrep(P::VPolygon{N})::AbstractHPolygon{N} where {N<:Real}
    error("this function is not implemented yet, see issue #5")
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
    return P.vertices_list
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
    if isempty(P.vertices_list)
        error("this polygon is empty")
    end
    i_max = 1
    @inbounds for i in 2:length(P.vertices_list)
        if dot(d, P.vertices_list[i] - P.vertices_list[i_max]) > zero(N)
            i_max = i
        end
    end
    return P.vertices_list[i_max]
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
    if isempty(P.vertices_list)
        error("this polygon is empty")
    end
    return P.vertices_list[1]
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
    if length(P.vertices_list) == 0
        return false
    elseif length(P.vertices_list) == 1
        return x == P.vertices_list[1]
    end

    zero_N = zero(N)
    if right_turn(P.vertices_list[1], x, P.vertices_list[end]) < zero_N
        return false
    end
    for i in 2:length(P.vertices_list)
        if right_turn(P.vertices_list[i], x, P.vertices_list[i-1]) < zero_N
            return false
        end
    end
    return true
end
