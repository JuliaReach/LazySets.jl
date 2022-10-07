export Polygon

"""
    Polygon{N, VN<:AbstractVector{N}} <: AbstractPolygon{N}

Type that represents a polygon by its line segments.

### Fields

- `line_segments` -- the list of line segments

### Notes

This type assumes that all line segments are sorted in counter-clockwise fashion.

To ensure this property, the constructor of `Polygon` runs a non-convex-hull
algorithm on the line segments by default. This also removes redundant line segments.
If the line segments are known to be sorted, the flag `apply_convex_hull=false` can
be used to skip this preprocessing.

### Examples

```
"""
struct Polygon{N, VN<:AbstractVector{N}} <: AbstractPolygon{N}
    line_segments::LineSegment{N, VN<:AbstractVector{N}}

    # default constructor that applies a convex hull algorithm
    function Polygon(line_segments::LineSegment{N, VN<:AbstractVector{N}};
                      apply_convex_hull::Bool=true,
                      algorithm::String="monotone_chain"
                     ) where {N, VN<:AbstractVector{N}}
        if apply_convex_hull
            vertices = convex_hull(vertices, algorithm=algorithm)
        end
        return new{N, VN}(line_segments)
    end
end
