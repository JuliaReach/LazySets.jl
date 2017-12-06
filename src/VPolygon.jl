import Base.<=

export VPolygon,
       vertices_list,
       singleton_list

"""
    VPolygon{N<:Real} <: LazySet

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
struct VPolygon{N<:Real} <: LazySet
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

"""
    dim(P::VPolygon)

Return the dimension of a polygon in vertex representation.

### Input

- `P` -- polygon in vertex representation

### Output

The ambient dimension of the polygon.
"""
function dim(P::VPolygon)
    2
end

"""
    σ(d::AbstractVector{<:Real}, P::VPolygon{N})::Vector{N} where {N<:Real}

Return the support vector of a polygon in a given direction.

### Input

- `d` -- direction
- `P` -- polygon in vertex representation

### Output

The support vector in the given direction. If the direction has norm zero, the
first vertex is returned.

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
function σ(d::AbstractVector{<:Real}, P::VPolygon{N})::Vector{N} where {N<:Real}
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

"""
    singleton_list(P::VPolygon{N})::Vector{Singleton{N}} where {N<:Real}

Return the vertices of a convex polygon in vertex representation as a list of
singletons.

### Input

- `P` -- a polygon vertex representation

### Output

List containing a singleton for each vertex.
"""
function singleton_list(P::VPolygon{N})::Vector{Singleton{N}} where {N<:Real}
    return [Singleton(vi) for vi in P.vertices_list]
end
