"""
    vertices(∅::EmptySet)

Construct an iterator over the vertices of an empty set.

### Input

- `∅` -- empty set

### Output

The empty iterator, as the empty set does not contain any vertices.
"""
function vertices(∅::EmptySet)
    N = eltype(∅)
    return EmptyIterator{Vector{N}}()
end
