function an_element(P::Polygon)
    @assert !isempty(P.vertices) "the polygon has no vertices"
    return P.vertices[1]
end
