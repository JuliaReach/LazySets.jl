"""
    vertices_list(P::VPolytope)

Return the list of vertices of a polytope in vertex representation.

### Input

- `P` -- polytope in vertex representation

### Output

The list of vertices.
"""
function vertices_list(P::VPolytope)
    return P.vertices
end
