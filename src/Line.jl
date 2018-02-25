export Line

"""
    Line{N<:Real} <: AbstractConvexSet{N}

Type that represents a line in 2D of the form ``a⋅x = b`` (i.e., a special case
of a `Hyperplane`).

### Fields

- `a` -- normal direction
- `b` -- constraint

### Examples

The line ``y = -x + 1``:

```jldoctest
julia> Line([1., 1.], 1.)
LazySets.Line{Float64}([1.0, 1.0], 1.0)
```
"""
struct Line{N<:Real} <: AbstractConvexSet{N}
    a::Vector{N}
    b::N

    # default constructor with length constraint
    Line{N}(a::Vector{N}, b::N) where {N<:Real} =
        (length(a) != 2 ? throw(DimensionMismatch) : new{N}(a, b))
end

# type-less convenience constructor
Line(a::Vector{N}, b::N) where {N<:Real} = Line{N}(a, b)

# constructor from a LinearConstraint
Line(c::LinearConstraint{N}) where {N<:Real} = Line{N}(c.a, c.b)


# --- LazySet interface functions ---


"""
    dim(L::Line)::Int

Return the ambient dimension of a line.

### Input

- `L` -- line

### Output

The ambient dimension of the line, which is 2.
"""
function dim(L::Line)::Int
    return 2
end

"""
    σ(d::AbstractVector{<:Real}, L::Line)::AbstractVector{<:Real}

Return the support vector of a line in a given direction.

### Input

- `d` -- direction
- `L` -- line

### Output

The support vector in the given direction, which is defined the same way as for
the more general `Hyperplane`.
"""
function σ(d::AbstractVector{<:Real}, L::Line)::AbstractVector{<:Real}
    return σ(d, Hyperplane(L.a, L.b))
end
