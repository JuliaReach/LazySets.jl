export VPolygonNC

"""
    VPolygonNC{N, VN<:AbstractVector{N}} <: LazySets{N}

Type that represents a non-convex polygon by its vertices.

### Fields

- `vertices` -- the list of vertices

### Notes

This type assumes that all vertices are sorted in counter-clockwise fashion.

To ensure this property, the constructor of `VPolygonNC` runs a convex-hull
algorithm on the vertices by default. This also removes redundant vertices.
If the vertices are known to be sorted, the flag `apply_convex_hull=false` can
be used to skip this preprocessing.

### Examples

```
"""
struct VPolygonNC{N, VN<:AbstractVector{N}} <: LazySets{N}
    edges::Vector{LineSegment{N, VN}}

    # default constructor that applies a convex hull algorithm
    function VPolygonNC(vertices::Vector{VN};
                        apply_convex_hull::Bool=true,
                        algorithm::String="monotone_chain"
                        ) where {N, VN<:AbstractVector{N}}
        if apply_convex_hull
            vertices = convex_hull(vertices, algorithm=algorithm)
        end
        return new{N, VN}(line_segments)
    end
end

# constructor with empty vertices list
VPolygonNC{N}() where {N} = VPolygon(Vector{Vector{N}}(), apply_convex_hull=false)
