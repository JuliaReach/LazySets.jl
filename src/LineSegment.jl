import Base.∈

export LineSegment

"""
    LineSegment{N<:Real} <: LazySet{N}

Type that represents a line segment in 2D between two points ``p`` and ``q``.

### Fields

- `p` -- first point
- `q` -- second point
"""
struct LineSegment{N<:Real} <: LazySet{N}
    p::AbstractVector{N}
    q::AbstractVector{N}

    # default constructor with length constraint
    LineSegment{N}(p::AbstractVector{N}, q::AbstractVector{N}) where {N<:Real} =
        (length(p) == length(q) == 2 ? new{N}(p, q) : throw(DimensionMismatch))
end

# type-less convenience constructor
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
    σ(d::AbstractVector{<:Real}, L::LineSegment)::AbstractVector{<:Real}

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
function σ(d::AbstractVector{<:Real}, L::LineSegment)::AbstractVector{<:Real}
    return sign(dot(L.q - L.p, d)) >= 0 ? L.q : L.p
end

"""
    ∈(x::AbstractVector{N}, L::LineSegment{N})::Bool where {N<:Real}

Check whether a given point is contained in a line segment.

### Input

- `x` -- point/vector
- `L` -- line segment

### Output

`true` iff ``x ∈ L``.

### Notes

The algorithm is inspired from [here](https://stackoverflow.com/a/328122).
"""
function ∈(x::AbstractVector{N}, L::LineSegment{N})::Bool where {N<:Real}
    @assert length(x) == dim(L)
    # check if the point is on the line through the line segment
    if abs((x[2] - L.p[2]) * (L.q[1] - L.p[1]) -
            (x[1] - L.p[1]) * (L.q[2] - L.p[2])) > 0
        return false
    end
    # check if the point is inside the box approximation of the line segment
    return min(L.p[1], L.q[1]) <= x[1] <= max(L.p[1], L.q[1]) &&
           min(L.p[2], L.q[2]) <= x[2] <= max(L.p[2], L.q[2])
end
