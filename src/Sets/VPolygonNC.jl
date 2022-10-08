export VPolygonNC

"""
    VPolygonNC{N, VN<:AbstractVector{N}} <: LazySets{N}

Type that represents a non-convex polygon by its vertices.

### Fields

- `vertices` -- the list of vertices

### Notes

This type assumes that all vertices are sorted in counter-clockwise fashion.

To ensure this property, the constructor of `VPolygonNC` applies an algorithm `initial_check` on the vertices by default. This also removes redundant vertices.


### Examples

```
"""
struct VPolygonNC{N, VN<:AbstractVector{N}} <: LazySets{N}
    edges::Vector{LineSegment{N, VN}}

    # default constructor that applies a convex hull algorithm
    function VPolygonNC(vertices::Vector{VN};
                        initial_check::Bool=true,
                        ) where {N, VN<:AbstractVector{N}}
        if initial_check
            vertices = convex_hull(vertices, algorithm=algorithm)
        end
        return new{N, VN}(line_segments)
    end
end

# constructor with empty vertices list
VPolygonNC{N}() where {N} = VPolygon(Vector{Vector{N}}(), apply_convex_hull=false)
