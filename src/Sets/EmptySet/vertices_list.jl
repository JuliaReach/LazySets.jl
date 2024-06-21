"""
    vertices_list(∅::EmptySet)

Return the list of vertices of an empty set.

### Input

- `∅` -- empty set

### Output

The empty list of vertices, as the empty set does not contain any vertices.
"""
function vertices_list(∅::EmptySet)
    N = eltype(∅)
    return Vector{Vector{N}}()
end
