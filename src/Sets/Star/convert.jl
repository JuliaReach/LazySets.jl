function convert(::Type{Star}, P::AbstractPolyhedron{N}) where {N}
    n = dim(P)
    c = zeros(N, n)
    V = Matrix(one(N) * I, n, n)
    return Star(c, V, P)
end
