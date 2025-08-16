"""
    vectorize(MZ::MatrixZonotope)

Vectorize a matrix zonotope to transform it into a zonotope.

### Input

- `MZ` -- a matrix zonotope

### Output

A zonotope.
"""
function vectorize(MZ::MatrixZonotope)
    c = vec(center(MZ))
    G = hcat([vec(Ai) for Ai in generators(MZ)]...)
    return Zonotope(c, G)
end

"""
    matrixize(Z::Zonotope, dims::Tuple{Int,Int})

Reshape a zonotope to transform it into a matrix zonotope.

### Input

- `Z`    -- a zonotope
- `dims` -- target dimensions of the matrix zonotope

### Output

A matrix zonotope.
"""
function matrixize(Z::Zonotope, dims::Tuple{Int,Int})
    @assert dims[1] * dims[2] == dim(Z) "cannot reshape a zonotope of dim = $(dim(Z)) into " *
                                        "a matrix zonotope of size $dims"

    c = Matrix(reshape(center(Z), dims))
    gens = [Matrix(reshape(col, dims)) for col in eachcol(genmat(Z))]
    return MatrixZonotope(c, gens)
end
