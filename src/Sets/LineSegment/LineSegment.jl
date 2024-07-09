"""
    LineSegment{N, VN<:AbstractVector{N}} <: AbstractZonotope{N}

Type that represents a line segment in 2D between two points ``p`` and ``q``.

### Fields

- `p` -- first point
- `q` -- second point

### Examples

A line segment along the ``x = y`` diagonal:

```jldoctest linesegment_constructor
julia> s = LineSegment([0., 0], [1., 1.])
LineSegment{Float64, Vector{Float64}}([0.0, 0.0], [1.0, 1.0])

julia> dim(s)
2
```

Use `plot(s)` to plot the extreme points of `s` and the line segment joining
them. If it is desired to remove the endpoints, pass the options
`markershape=:none` and `seriestype=:shape`.

Membership is checked with ∈ (`in`):

```jldoctest linesegment_constructor
julia> [0., 0] ∈ s && [.25, .25] ∈ s && [1., 1] ∈ s && [.5, .25] ∉ s
true
```

We can check whether the intersection with another line segment is empty, and
optionally compute a witness (which is the unique common point in this case):

```jldoctest linesegment_constructor
julia> sn = LineSegment([1., 0], [0., 1.])
LineSegment{Float64, Vector{Float64}}([1.0, 0.0], [0.0, 1.0])

julia> isdisjoint(s, sn)
false

julia> isdisjoint(s, sn, true)
(false, [0.5, 0.5])
```
"""
struct LineSegment{N,VN<:AbstractVector{N}} <: AbstractZonotope{N}
    p::VN
    q::VN

    # default constructor with length constraint
    function LineSegment(p::VN, q::VN) where {N,VN<:AbstractVector{N}}
        @assert length(p) == length(q) == 2 "points for line segments must " *
                                            "be two-dimensional, but their lengths are $(length(p)) and " *
                                            "$(length(q))"
        return new{N,VN}(p, q)
    end
end
