"""
    MatrixZonotopeExp{N}(MZ::AbstractMatrixZonotope{N}) <: AbstractMatrixZonotope{N}

Represents the matrix exponential of a matrix zonotope.

### Fields

- `M` -- a square matrix zonotope representing the exponent.

### Notes

Mathematically, this represents the set:

```math
\\mathcal{C} = \\{ \\exp(M) ~|~ M \\in \\mathcal{M} \\}
```
"""
struct MatrixZonotopeExp{N,T<:AbstractMatrixZonotope{N}} <: AbstractMatrixZonotope{N}
    M::T

    function MatrixZonotopeExp(MZ::T) where {N,T<:AbstractMatrixZonotope{N}}
        @assert size(MZ, 1) == size(MZ, 2) "the exponent must be a square matrix zonotope"
        return new{N,T}(MZ)
    end
end

Base.size(EMZ::MatrixZonotopeExp) = size(EMZ.M)
Base.size(EMZ::MatrixZonotopeExp, d::Int) = size(EMZ.M, d)
