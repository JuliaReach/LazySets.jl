export LinearConstraint, Line, intersection

"""
    LinearConstraint

Type that represents a linear constraint (a half-space) of the form a⋅x ≦ b.

### Fields

- `a` --  a normal direction
- `b` -- the constraint

EXAMPLES:

The set `y >= 0` in the plane:

```julia
julia> LinearConstraint([0, -1.], 0.)
LazySets.LinearConstraint([0.0, -1.0], 0.0)
```
"""
struct LinearConstraint{N<:Real}
    a::Vector{N}
    b::N
end

"""
    Line

Type that represents a line in 2D of the form a⋅x = b.

### Fields

- `a`  -- a normal direction (size = 2)
- `b`  -- the constraint

### Examples

The line `y = -x + 1`:

```julia
julia> Line([1., 1.], 1.)
LazySets.Line([1.0, 1.0], 1.0)
```
"""
struct Line{N<:Real}
    a::Vector{N}
    b::N

    # default constructor with length constraint
    Line{N}(a::Vector{N}, b::N) where {N<:Real} =
        (length(a) != 2
            ? throw(DimensionMismatch)
            : new{N}(a, b))
end
# type-less convenience constructor
Line(a::Vector{N}, b::N) where {N<:Real} = Line{N}(a, b)
# constructor from a LinearConstraint
Line(c::LinearConstraint{N}) where {N<:Real} = Line{N}(c.a, c.b)

"""
    intersection(Δ1, Δ2)

Return the intersection of two 2D lines.

### Input

- `Δ1` -- a line

- `Δ2` -- another line

### Output

The intersection point.

### Examples

The line `y = -x + 1` intersected with `y = x`:

```julia
julia> intersection(Line([1., 1.], 1.), Line([-1., 1.], 0.))
2-element Array{Float64,1}:
 0.5
 0.5
```
"""
function intersection(Δ1::Line{N}, Δ2::Line{N})::Vector{N} where {N<:Real}
    b = [Δ1.b, Δ2.b]
    a = [Δ1.a.' ; Δ2.a.']
    return a \ b
end

