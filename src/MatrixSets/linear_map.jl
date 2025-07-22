"""
    linear_map(M::AbstractMatrix, MZ::MatrixZonotope)

Applies a linear transformation to a matrix zonotope from the left.

### Input
- `M` -- a linear map / matrix
- `MZ`` -- a matrix zonotope

### Output

A matrix zonotope with transformed center and generators
"""
function linear_map(M::AbstractMatrix, MZ::MatrixZonotope)
    @assert size(M, 2) == size(MZ, 1) "incompatible dimensions"
    gens = [M * Ai for Ai in generators(MZ)]
    return MatrixZonotope(M * center(MZ), gens)
end

"""
    linear_map(M::AbstractMatrix, MZ::MatrixZonotope)

Applies a linear transformation to a matrix zonotope from the right.

### Input
- `M` -- a linear map / matrix
- `MZ`` -- a matrix zonotope

### Output

A matrix zonotope with transformed center and generators
"""
function linear_map(MZ::MatrixZonotope, M::AbstractMatrix)
    @assert size(MZ, 2) == size(M, 1) "incompatible dimensions"
    gens = [Ai * M for Ai in generators(MZ)]
    return MatrixZonotope(center(MZ) * M, gens)
end
