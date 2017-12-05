export LinearConstraint,
       Line,
       intersection

"""
    LinearConstraint{N<:Real}

Type that represents a linear constraint (a half-space) of the form ``a⋅x ≤ b``.

### Fields

- `a` -- normal direction
- `b` -- constraint

### Examples

The set ``y ≥ 0`` in the plane:

```jldoctest
julia> LinearConstraint([0, -1.], 0.)
LazySets.LinearConstraint{Float64}([0.0, -1.0], 0.0)
```
"""
struct LinearConstraint{N<:Real}
    a::Vector{N}
    b::N
end

"""
    Line{N<:Real}

Type that represents a line in 2D of the form ``a⋅x = b``.

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
struct Line{N<:Real}
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

"""
    intersection(L1::Line{N}, L2::Line{N})::Vector{N} where {N<:Real}

Return the intersection of two 2D lines.

### Input

- `L1` -- first line

- `L2` -- second line

### Output

The intersection point, if any.
Throws a `SingularException` if the lines do not intersect.

### Examples

The line ``y = -x + 1`` intersected with the line ``y = x``:

```jldoctest
julia> intersection(Line([-1., 1.], 0.), Line([1., 1.], 1.))
2-element Array{Float64,1}:
 0.5
 0.5
```
"""
function intersection(L1::Line{N}, L2::Line{N})::Vector{N} where {N<:Real}
    b = [L1.b, L2.b]
    a = [L1.a.'; L2.a.']
    return a \ b
end
