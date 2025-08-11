"""
    vectorize(MZ::MatrixZonotope)

Vectorize a matrix zonotope to transform it into a zonotope.

### Input

- `MZ`       -- a matrix zonotope

### Output

A zonotope.
"""
function vectorize(MZ::MatrixZonotope)
    c = vec(center(MZ))
    G = hcat([vec(Ai) for Ai in generators(MZ)]...)
    return Zonotope(c, G)
end

"""
    matrixize(Z::Zonotope)

Reshape a zonotope zonotope to transform it into a matrix zonotope.

### Input

- `Z`       -- a zonotope

### Output

A matrix zonotope.
"""
function matrixize(Z::Zonotope, dim::Tuple{Int,Int})
    @assert dim[1] * dim[2] == dim(Z) "cannot reshape a zonotope of dim = $(dim(Z)) into " *
                                      "a matrix zonotope of size $dim"

    c = reshape(center(Z), dim)
    gens = [Matrix(reshape(col, dim)) for col in eachcol(genmat(G))]
    return MatrixZonotope(c, gens)
end
