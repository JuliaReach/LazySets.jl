"""
    AbstractMatrixZonotope{N}

Abstract supertype for all matrix zonotope representations.

Every concrete `AbstractSparsePolynomialZonotope` must define the following functions:
- `size(::AbstractSparsePolynomialZonotope)` -- return the size of the matrix zonotope
"""
abstract type AbstractMatrixZonotope{N} end

"""
    size(::AbstractMatrixZonotope, [dim])

Return a tuple containing the dimensions of a matrix zonotope. Optionally you can specify a dimension to just get the length of that dimension.
"""
function Base.size(::AbstractMatrixZonotope) end
