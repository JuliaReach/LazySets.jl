"""
# Extended help

    translate(L::LineSegment, v::AbstractVector)

### Algorithm

We add the vector to both defining points of the line segment.
"""
@validate function translate(L::LineSegment, v::AbstractVector)
    return translate!(deepcopy(L), v)
end

@validate function translate!(L::LineSegment, v::AbstractVector)
    L.p .+= v
    L.q .+= v
    return L
end
