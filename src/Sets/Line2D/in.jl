"""
    ∈(x::AbstractVector, L::Line2D)

Check whether a given point is contained in a 2D line.

### Input

- `x` -- point/vector
- `L` -- 2D line

### Output

`true` iff `x ∈ L`.

### Algorithm

The point ``x`` belongs to the line if and only if ``a⋅x = b`` holds.
"""
function ∈(x::AbstractVector, L::Line2D)
    @assert length(x) == 2 "a $(length(x))-dimensional vector is " *
                           "incompatible with a 2-dimensional line"
    return _isapprox(dot(L.a, x), L.b)
end
