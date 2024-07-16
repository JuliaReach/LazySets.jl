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
