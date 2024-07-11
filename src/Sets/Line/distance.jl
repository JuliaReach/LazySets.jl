"""
    distance(x::AbstractVector, L::Line; [p]::Real=2.0)

Compute the distance between point `x` and the line with respect to the given
`p`-norm.

### Input

- `x` -- point/vector
- `L` -- line
- `p` -- (optional, default: `2.0`) the `p`-norm used; `p = 2.0` corresponds to
         the usual Euclidean norm

### Output

A scalar representing the distance between `x` and the line `L`.
"""
@commutative function distance(x::AbstractVector, L::Line; p::Real=2.0)
    d = L.d  # direction of the line
    t = dot(x - L.p, d) / dot(d, d)
    return distance(x, L.p + t * d; p=p)
end
