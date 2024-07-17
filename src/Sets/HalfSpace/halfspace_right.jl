"""
    halfspace_right(p::AbstractVector, q::AbstractVector)

Return a half-space describing the 'right' of a two-dimensional line segment
through two points.

### Input

- `p` -- first point
- `q` -- second point

### Output

The half-space whose boundary goes through the two points `p` and `q` and which
describes the right-hand side of the directed line segment `pq`.

### Algorithm

See the documentation of [`halfspace_left`](@ref).
"""
function halfspace_right(p::AbstractVector, q::AbstractVector)
    return halfspace_left(q, p)
end
