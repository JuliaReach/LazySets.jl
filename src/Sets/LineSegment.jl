export LineSegment,
       halfspace_left, halfspace_right

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

isoperationtype(::Type{<:LineSegment}) = false

"""
    dim(L::LineSegment)

Return the ambient dimension of a 2D line segment.

### Input

- `L` -- 2D line segment

### Output

The ambient dimension of the 2D line segment, which is ``2``.
"""
function dim(L::LineSegment)
    return 2
end

"""
    σ(d::AbstractVector, L::LineSegment)

Return the support vector of a 2D line segment in a given direction.

### Input

- `d` -- direction
- `L` -- 2D line segment

### Output

The support vector in the given direction.

### Algorithm

If the angle between the vector ``q - p`` and ``d`` is bigger than 90° and less
than 270° (measured in counter-clockwise order), the result is ``p``, otherwise
it is ``q``.
If the angle is exactly 90° or 270°, or if the direction has norm zero, this
implementation returns ``q``.
"""
function σ(d::AbstractVector, L::LineSegment)
    return sign(dot(L.q - L.p, d)) >= 0 ? L.q : L.p
end

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

"""
    an_element(L::LineSegment)

Return some element of a 2D line segment.

### Input

- `L` -- 2D line segment

### Output

The first vertex of the line segment.
"""
function an_element(L::LineSegment)
    return L.p
end

"""
    ∈(x::AbstractVector, L::LineSegment)

Check whether a given point is contained in a 2D line segment.

### Input

- `x` -- point/vector
- `L` -- 2D line segment

### Output

`true` iff ``x ∈ L``.

### Algorithm

Let ``L = (p, q)`` be the line segment with extreme points ``p`` and ``q``, and
let ``x`` be the given point.

1. A necessary condition for ``x ∈ (p, q)`` is that the three points are
   aligned, thus their cross product should be zero.
2. It remains to check that ``x`` belongs to the box approximation of ``L``.
   This amounts to comparing each coordinate with those of the extremes ``p``
   and ``q``.

### Notes

The algorithm is inspired from [here](https://stackoverflow.com/a/328110).
"""
function ∈(x::AbstractVector, L::LineSegment)
    @assert length(x) == dim(L) "a vector of length $(length(x)) is " *
                                "incompatible with a $(dim(L))-dimensional set"

    # check if point x is on the line through the line segment (p, q)
    p = L.p
    q = L.q
    if isapproxzero(right_turn(p, q, x))
        # check if the point is inside the box approximation of the line segment
        return @inbounds (_leq(min(p[1], q[1]), x[1]) &&
                          _leq(x[1], max(p[1], q[1])) &&
                          _leq(min(p[2], q[2]), x[2]) &&
                          _leq(x[2], max(p[2], q[2])))
    else
        return false
    end
end

"""
    center(L::LineSegment)

Return the center of a 2D line segment.

### Input

- `L` -- 2D line segment

### Output

The center of the line segment.
"""
function center(L::LineSegment)
    return L.p + (L.q - L.p) / 2
end

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

"""
    vertices_list(L::LineSegment)

Return the list of vertices of a 2D line segment.

### Input

- `L` -- 2D line segment

### Output

The list of end points of the line segment.
"""
function vertices_list(L::LineSegment)
    return [L.p, L.q]
end

"""
    rand(::Type{LineSegment}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

Create a random 2D line segment.

### Input

- `LineSegment` -- type for dispatch
- `N`           -- (optional, default: `Float64`) numeric type
- `dim`         -- (optional, default: 2) dimension
- `rng`         -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`        -- (optional, default: `nothing`) seed for reseeding

### Output

A random 2D line segment.

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.
"""
function rand(::Type{LineSegment};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing)
    @assert dim == 2 "cannot create a random LineSegment of dimension $dim"
    rng = reseed!(rng, seed)
    p = randn(rng, N, dim)
    q = randn(rng, N, dim)
    return LineSegment(p, q)
end

"""
    halfspace_left(L::LineSegment)

Return a half-space describing the 'left' of a two-dimensional 2D line segment
through two points.

### Input

- `L` -- 2D line segment

### Output

The half-space whose boundary goes through the two points `p` and `q` and which
describes the left-hand side of the directed line segment `pq`.
"""
halfspace_left(L::LineSegment) = halfspace_left(L.p, L.q)

"""
    halfspace_right(L::LineSegment)

Return a half-space describing the 'right' of a two-dimensional 2D line segment
through two points.

### Input

- `L` -- 2D line segment

### Output

The half-space whose boundary goes through the two points `p` and `q` and which
describes the right-hand side of the directed line segment `pq`.
"""
halfspace_right(L::LineSegment) = halfspace_right(L.p, L.q)

"""
    constraints_list(L::LineSegment)

Return a list of constraints defining a 2D line segment in 2D.

### Input

- `L` -- 2D line segment

### Output

A vector of constraints defining the line segment.

### Algorithm

``L`` is defined by 4 constraints. In this algorithm, the first two constraints
are returned by ``halfspace_right`` and ``halfspace_left``, and the other two
are obtained by considering a vector parallel to the line segment passing
through one of the vertices.
"""
function constraints_list(L::LineSegment)
    p, q = L.p, L.q
    d = @inbounds [p[2] - q[2], q[1] - p[1]]
    return [halfspace_left(L), halfspace_right(L),
            halfspace_right(p, p + d), halfspace_left(q, q + d)]
end

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

"""
    ngens(L::LineSegment)

Return the number of generators of a 2D line segment.

### Input

- `L` -- 2D line segment

### Output

The number of generators.

### Algorithm

A line segment has either one generator, or zero generators if it is a
degenerated line segment of length zero.
"""
function ngens(L::LineSegment)
    return _isapprox(L.p, L.q) ? 0 : 1
end

function scale!(α::Real, L::LineSegment)
    L.p .*= α
    L.q .*= α
    return L
end
