"""
    halfspace_left(p::AbstractVector, q::AbstractVector)

Return a half-space describing the 'left' of a two-dimensional line segment
through two points.

### Input

- `p` -- first point
- `q` -- second point

### Output

The half-space whose boundary goes through the two points `p` and `q` and which
describes the left-hand side of the directed line segment `pq`.

### Algorithm

The half-space ``a⋅x ≤ b`` is calculated as `a = [dy, -dx]`, where
``d = (dx, dy)`` denotes the line segment `pq`, i.e.,
``\\vec{d} = \\vec{p} - \\vec{q}``, and `b = dot(p, a)`.

### Examples

The left half-space of the "east" and "west" directions in two-dimensions are
the upper and lower half-spaces:

```jldoctest halfspace_left
julia> using LazySets: halfspace_left

julia> halfspace_left([0.0, 0.0], [1.0, 0.0])
HalfSpace{Float64, Vector{Float64}}([0.0, -1.0], 0.0)

julia> halfspace_left([0.0, 0.0], [-1.0, 0.0])
HalfSpace{Float64, Vector{Float64}}([0.0, 1.0], 0.0)
```

We create a box from the sequence of line segments that describe its edges:

```jldoctest halfspace_left
julia> H1 = halfspace_left([-1.0, -1.0], [1.0, -1.0]);

julia> H2 = halfspace_left([1.0, -1.0], [1.0, 1.0]);

julia> H3 = halfspace_left([1.0, 1.0], [-1.0, 1.0]);

julia> H4 = halfspace_left([-1.0, 1.0], [-1.0, -1.0]);

julia> H = HPolygon([H1, H2, H3, H4]);

julia> B = BallInf([0.0, 0.0], 1.0);

julia> B ⊆ H && H ⊆ B
true
```
"""
function halfspace_left(p::AbstractVector, q::AbstractVector)
    @assert length(p) == length(q) == 2 "the points must be two-dimensional"
    @assert p != q "the points must not be equal"
    a = [q[2] - p[2], p[1] - q[1]]
    return HalfSpace(a, dot(p, a))
end
