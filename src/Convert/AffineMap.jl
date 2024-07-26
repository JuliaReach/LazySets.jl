"""
    convert(::Type{STAR}, P::AbstractPolyhedron{N}) where {N}

Convert a polyhedral set to a star set represented as a lazy affine map.

### Input

- `STAR` -- target type
- `P`    -- polyhedral set

### Output

A star set.
"""
function convert(::Type{STAR}, P::AbstractPolyhedron{N}) where {N}
    n = dim(P)
    c = zeros(N, n)
    V = Matrix(one(N) * I, n, n)
    return AffineMap(V, P, c)
end

"""
    convert(::Type{STAR}, X::Star)

Convert a star set to its equivalent representation as a lazy affine map.

### Input

- `STAR` -- target type
- `X`    -- star set

### Output

A star set.
"""
function Base.convert(::Type{STAR}, X::Star)
    return AffineMap(X.V, X.P, X.c)
end
