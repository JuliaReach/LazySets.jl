"""
   genmat(L::LineSegment)

Return the generator matrix of a 2D line segment.

### Input

- `L` -- 2D line segment

### Output

A matrix with at most one column representing the generator of `L`.
"""
function genmat(L::LineSegment)
    N = eltype(L)
    if _isapprox(L.p, L.q)
        # degenerate line segment has no generators
        return zeros(N, dim(L), 0)
    end
    return hcat((L.p - L.q) / 2)
end
