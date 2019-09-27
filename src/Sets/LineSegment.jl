import Base: rand,
             ∈

export LineSegment,
       halfspace_left, halfspace_right,
       constraints_list

"""
    LineSegment{N<:Real} <: AbstractZonotope{N}

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

Use `plot(s)` to plot the extreme points of `s` and the line segment joining
them. Membership test is computed with ∈ (`in`):

```jldoctest linesegment_constructor
julia> [0., 0] ∈ s && [.25, .25] ∈ s && [1., 1] ∈ s && [.5, .25] ∉ s
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
struct LineSegment{N<:Real} <: AbstractZonotope{N}
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

isoperationtype(::Type{<:LineSegment}) = false

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


# --- AbstractZonotope interface functions ---


"""
   genmat(L::LineSegment)

Return the generator matrix of a line segment.

### Input

- `L` -- line segment

### Output

A matrix with a single column representing the generator of `L`.
"""
function genmat(L::LineSegment)
    return genmat_fallback(L)
end

"""
    generators(L::LineSegment{N}) where {N<:Real}

Return an iterator over the (single) generator of a line segment.

### Input

- `L` -- line segment

### Output

A one-element iterator over the generator of `L`.
"""
function generators(L::LineSegment{N}) where {N<:Real}
    if L.p == L.q
        # degenerate line segment has no generators
        return EmptyGeneratorIterator{N}()
    end
    return [L.p - center(L)]
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


# --- LazySet interface functions ---


"""
    rand(::Type{LineSegment}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing
        )::LineSegment{N}

Create a random line segment.

### Input

- `LineSegment` -- type for dispatch
- `N`           -- (optional, default: `Float64`) numeric type
- `dim`         -- (optional, default: 2) dimension
- `rng`         -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`        -- (optional, default: `nothing`) seed for reseeding

### Output

A random line segment.

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.
"""
function rand(::Type{LineSegment};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int, Nothing}=nothing
             )::LineSegment{N}
    @assert dim == 2 "cannot create a random LineSegment of dimension $dim"
    rng = reseed(rng, seed)
    p = randn(rng, N, dim)
    q = randn(rng, N, dim)
    return LineSegment(p, q)
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

"""
    constraints_list(L::LineSegment{N}) where {N<:Real}

Return the list of constraints defining a line segment in 2D.

### Input

- `L` -- line segment

### Output

A vector of constraints that define the line segment.

### Algorithm

``L`` is defined by 4 constraints. In this algorithm, the first two constraints
are returned by ``halfspace_right`` and ``halfspace_left``, and the other two
are obtained by considering the vector normal to the line segment that passes
through each opposite vertex.

### Notes

This function returns a vector of halfspaces. It does not return equality
constraints.
"""
function constraints_list(L::LineSegment{N}) where {N<:Real}
    p, q = L.p, L.q
    d = [(p[2]-q[2]), (q[1]-p[1])]
    return [halfspace_left(L), halfspace_right(L),
            halfspace_right(p, p + d), halfspace_left(q, q + d)]
end

"""
    translate(L::LineSegment{N}, v::AbstractVector{N}) where {N<:Real}

Translate (i.e., shift) a line segment by a given vector.

### Input

- `L` -- line segment
- `v` -- translation vector

### Output

A translated line segment.

### Algorithm

We add the vector to both defining points of the line segment.
"""
function translate(L::LineSegment{N}, v::AbstractVector{N}) where {N<:Real}
    @assert length(v) == dim(L) "cannot translate a $(dim(L))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    return LineSegment(L.p + v, L.q + v)
end

"""
    plot_recipe(L::LineSegment{N}, [ε]::N=zero(N)) where {N<:Real}

Convert a line segment to a pair `(x, y)` of points for plotting.

### Input

- `L` -- line segment
- `ε` -- (optional, default: `0`) ignored, used for dispatch

### Output

A pair `(x, y)` of two points that can be plotted.
"""
function plot_recipe(L::LineSegment{N}, ε::N=zero(N)) where {N<:Real}
    return [L.p[1], L.q[1]], [L.p[2], L.q[2]]
end
