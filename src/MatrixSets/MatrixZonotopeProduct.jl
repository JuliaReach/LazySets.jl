"""
    MatrixZonotopeProduct{N, MAT1<:MatrixZonotope{N}, NM,
        MAT2<:MatrixZonotope{NM}}(A::MAT1, B::MAT2) <: AbstractMatrixZonotope{N}

Represents the product of two matrix zonotopes.

### Fields

- `A` -- a matrix zonotope representing the left factor in the product
- `B` -- a matrix zonotope representing the right factor in the product

### Notes

Mathematically, this represents the set of matrices:

```math
\\mathcal{C} = \\{ A \\cdot B ~|~ A \\in \\mathcal{A},\\ B \\in \\mathcal{B} \\}
```
"""
struct MatrixZonotopeProduct{N, MAT1<:MatrixZonotope{N}, NM, MAT2<:MatrixZonotope{NM}} <: AbstractMatrixZonotope{N}
    A::MAT1
    B::MAT2

    function MatrixZonotopeProduct(A::MAT1, B::MAT2) where {N, MAT1<:MatrixZonotope{N}, NM, MAT2<:MatrixZonotope{NM}}
        @assert size(A.A0, 2) == size(B.A0, 1) "incompatible dimensions"
        return new{N, MAT1, NM, MAT2}(A, B)
    end
end

"""
    *(A::MatrixZonotope, B::MatrixZonotope)

Alias to create a `MatrixZonotopeProduct` object.
"""
function Base.:*(A::MatrixZonotope, B::MatrixZonotope)
    return MatrixZonotopeProduct(A, B)
end

Base.size(M::MatrixZonotopeProduct) = (size(M.A, 1), size(M.B, 2))
Base.size(M::MatrixZonotopeProduct, d::Int) = size(M)[d]
