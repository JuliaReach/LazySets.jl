"""
    convex_hull(P::VPolygon, Q::VPolygon; [algorithm]::String="monotone_chain")

Return the convex hull of two polygons in vertex representation.

### Input

- `P`         -- polygon in vertex representation
- `Q`         -- polygon in vertex representation
- `algorithm` -- (optional, default: "monotone_chain") the algorithm used to
                 compute the convex hull

### Output

A new polygon such that its vertices are the convex hull of the two polygons.

### Notes

The vertices of the output polygon are sorted in counter-clockwise fashion.
"""
function convex_hull(P::VPolygon, Q::VPolygon;
                     algorithm::String="monotone_chain")
    vunion = [P.vertices; Q.vertices]
    convex_hull!(vunion; algorithm=algorithm)
    return VPolygon(vunion; apply_convex_hull=false)
end
