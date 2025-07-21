"""
    AbstractMatrixZonotope{N}

Abstract supertype for all matrix zonotope representations, including standard (`MatrixZonotope`), products (`MatrixZonotopeProduct`), polynomial matrix zonotopes, and others.

Every concrete `AbstractSparsePolynomialZonotope` must define the following functions:
- `size(::AbstractSparsePolynomialZonotope)` -- return the size of the matrix zonotope
"""
abstract type AbstractMatrixZonotope{N} end

"""
    size(::AbstractMatrixZonotope)

Return the dimensions of the matrices in the matrix zonotope set. 

"""
function Base.size(::AbstractMatrixZonotope) end

"""
    size(::AbstractMatrixZonotope, ::Int)

Return the dimensions of a matrix zonotope along a specified axis.  

"""
function Base.size(::AbstractMatrixZonotope, ::Int) end
