import Base.∈

export LineSegment,
       halfspace_left, halfspace_right

"""
    LineSegment{N<:Real} <: AbstractCentrallySymmetricPolytope{N}

Type that represents a line segment in 2D between two points ``p`` and ``q``.

### Fields

- `p` -- first point
- `q` -- second point

### Examples

A line segment along the ``x = y`` diagonal:

```jldoctest linesegment_constructor
julia> s = LineSegment([0., 0], [1., 1.])
LineSegment{Float64}([0.0, 0.0], [1.0, 1.0])
julia> dim(s)
2
```

Use ``plot(s)`` to plot the extreme points of `s` and the line segment joining
them.
Membership test is computed with ∈ (`in`):

```jldoctest linesegment_constructor
julia> [0., 0] ∈ s && [.25, .25] ∈ s && [1., 1] ∈ s && !([.5, .25] ∈ s)
true
```

We can check the intersection with another line segment, and optionally compute
a witness (which is just the common point in this case):

```jldoctest linesegment_constructor
julia> sn = LineSegment([1., 0], [0., 1.])
LineSegment{Float64}([1.0, 0.0], [0.0, 1.0])
julia> isempty(s ∩ sn)
false
julia> is_intersection_empty(s, sn, true)
(false, [0.5, 0.5])
```
"""
struct LineSegment{N<:Real} <: AbstractCentrallySymmetricPolytope{N}
    p::AbstractVector{N}
    q::AbstractVector{N}

    # default constructor with length constraint
    function LineSegment{N}(p::AbstractVector{N},
                            q::AbstractVector{N}) where {N<:Real}
        @assert length(p) == length(q) == 2 "points for line segments must " *
            "be two-dimensional"
        return new{N}(p, q)
    end
end

# convenience constructor without type parameter
LineSegment(p::AbstractVector{N}, q::AbstractVector{N}) where {N<:Real} =
    LineSegment{N}(p, q)


# --- LazySet interface functions ---


"""
    dim(L::LineSegment)::Int

Return the ambient dimension of a line segment.

### Input

- `L` -- line segment

### Output

The ambient dimension of the line segment, which is 2.
"""
function dim(L::LineSegment)::Int
    return 2
end

"""
    σ(d::AbstractVector{N}, L::LineSegment{N}) where {N<:Real}

Return the support vector of a line segment in a given direction.

### Input

- `d` -- direction
- `L` -- line segment

### Output

The support vector in the given direction.

### Algorithm

If the angle between the vector ``q - p`` and ``d`` is bigger than 90° and less
than 270° (measured in counter-clockwise order), the result is ``p``, otherwise
it is ``q``.
If the angle is exactly 90° or 270°, or if the direction has norm zero, this
implementation returns ``q``.
"""
function σ(d::AbstractVector{N}, L::LineSegment{N}) where {N<:Real}
    return sign(dot(L.q - L.p, d)) >= 0 ? L.q : L.p
end

"""
    an_element(L::LineSegment{N}) where {N<:Real}

Return some element of a line segment.

### Input

- `L` -- line segment

### Output

The first vertex of the line segment.
"""
function an_element(L::LineSegment{N}) where {N<:Real}
    return L.p
end

"""
    ∈(x::AbstractVector{N}, L::LineSegment{N})::Bool where {N<:Real}

Check whether a given point is contained in a line segment.

### Input

- `x` -- point/vector
- `L` -- line segment

### Output

`true` iff ``x ∈ L``.

### Algorithm

Let ``L = (p, q)`` be the line segment with extremes ``p`` and ``q``, and let
``x`` be the given point.

1. A necessary condition for ``x ∈ (p, q)`` is that the three points are aligned,
   thus their cross product should be zero.
2. It remains to check that ``x`` belongs to the box approximation of ``L``.
   This amounts to comparing each coordinate with those of the extremes ``p``
   and ``q``.

### Notes

The algorithm is inspired from [here](https://stackoverflow.com/a/328110).
"""
function ∈(x::AbstractVector{N}, L::LineSegment{N})::Bool where {N<:Real}
    @assert length(x) == dim(L)
    # check if the point is on the line through the line segment
    if (x[2] - L.p[2]) * (L.q[1] - L.p[1]) -
            (x[1] - L.p[1]) * (L.q[2] - L.p[2]) != 0
        return false
    end
    # check if the point is inside the box approximation of the line segment
    return min(L.p[1], L.q[1]) <= x[1] <= max(L.p[1], L.q[1]) &&
           min(L.p[2], L.q[2]) <= x[2] <= max(L.p[2], L.q[2])
end


# --- AbstractCentrallySymmetric interface functions ---


"""
    center(L::LineSegment{N})::Vector{N} where {N<:Real}

Return the center of a line segment.

### Input

- `L` -- line segment

### Output

The center of the line segment.
"""
function center(L::LineSegment{N})::Vector{N} where {N<:Real}
    return L.p + (L.q - L.p) / 2
end


# --- AbstractPolytope interface functions ---


"""
    vertices_list(L::LineSegment{N}
                 )::Vector{<:AbstractVector{N}} where {N<:Real}

Return the list of vertices of a line segment.

### Input

- `L` -- line segment

### Output

The list of end points of the line segment.
"""
function vertices_list(L::LineSegment{N}
                      )::Vector{<:AbstractVector{N}} where {N<:Real}
    return [L.p, L.q]
end


# --- LineSegment functions ---


"""
    halfspace_left(L::LineSegment)

Return a half-space describing the 'left' of a two-dimensional line segment
through two points.

### Input

 - `L` -- line segment

### Output

The half-space whose boundary goes through the two points `p` and `q` and which
describes the left-hand side of the directed line segment `pq`.
"""
halfspace_left(L::LineSegment) = halfspace_left(L.p, L.q)

"""
    halfspace_right(L::LineSegment)

Return a half-space describing the 'right' of a two-dimensional line segment
through two points.

### Input

 - `L` -- line segment

### Output

The half-space whose boundary goes through the two points `p` and `q` and which
describes the right-hand side of the directed line segment `pq`.
"""
halfspace_right(L::LineSegment) = halfspace_right(L.p, L.q)
