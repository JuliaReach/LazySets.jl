"""
    translate(L::LineSegment, v::AbstractVector)

Translate (i.e., shift) a 2D line segment by a given vector.

### Input

- `L` -- 2D line segment
- `v` -- translation vector

### Output

A translated line segment.

### Algorithm

We add the vector to both defining points of the line segment.
"""
function translate(L::LineSegment, v::AbstractVector)
    @assert length(v) == dim(L) "cannot translate a $(dim(L))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    return LineSegment(L.p + v, L.q + v)
end
