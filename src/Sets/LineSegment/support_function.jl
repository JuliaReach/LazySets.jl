"""
    ρ(d::AbstractVector, L::LineSegment)

Evaluate the support function of a 2D line segment in a given direction.

### Input

- `d` -- direction
- `L` -- 2D line segment

### Output

Evaluation of the support function in the given direction.
"""
function ρ(d::AbstractVector, L::LineSegment)
    return max(dot(L.p, d), dot(L.q, d))
end
