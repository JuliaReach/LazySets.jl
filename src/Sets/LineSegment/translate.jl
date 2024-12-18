"""
# Extended help

    translate(L::LineSegment, v::AbstractVector)

### Algorithm

We add the vector to both defining points of the line segment.
"""
function translate(L::LineSegment, v::AbstractVector)
    @assert length(v) == dim(L) "cannot translate a $(dim(L))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    return LineSegment(L.p + v, L.q + v)
end
