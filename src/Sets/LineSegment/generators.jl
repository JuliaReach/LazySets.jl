"""
    generators(L::LineSegment)

Return an iterator over the (single) generator of a 2D line segment.

### Input

- `L` -- 2D line segment

### Output

An iterator over the generator of `L`, if any.
"""
function generators(L::LineSegment)
    if _isapprox(L.p, L.q)
        # degenerate line segment has no generators
        N = eltype(L)
        return EmptyIterator{Vector{N}}()
    end
    return SingletonIterator((L.p - L.q) / 2)
end
