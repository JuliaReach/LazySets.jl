"""
    convert(::Type{Star}, P::AbstractPolyhedron{N}) where {N}

Convert a polyhedral set to a star set.

### Input

- `Star` -- target type
- `P`    -- polyhedral set

### Output

A star set.
"""
function convert(::Type{Star}, P::AbstractPolyhedron{N}) where {N}
    n = dim(P)
    c = zeros(N, n)
    V = Matrix(one(N) * I, n, n)
    return Star(c, V, P)
end
