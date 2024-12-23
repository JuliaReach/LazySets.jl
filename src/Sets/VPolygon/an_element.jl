"""
# Extended help

    an_element(P::VPolygon)

### Output

The first vertex of the polygon in vertex representation.
"""
function an_element(P::VPolygon)
    @assert !isempty(P.vertices) "the polygon has no vertices"
    return P.vertices[1]
end
