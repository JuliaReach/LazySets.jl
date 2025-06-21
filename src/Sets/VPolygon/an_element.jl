function an_element(P::VPolygon)
    @assert !isempty(P.vertices) "the polygon has no vertices"
    return P.vertices[1]
end
