export inner

"""
    inner(x::AbstractVector{N}, A::AbstractMatrix{N}, y::AbstractVector{N}
         ) where {N}

Compute the inner product ``xᵀ A y``.

### Input

- `x` -- vector on the left
- `A` -- matrix
- `y` -- vector on the right

### Output

The (scalar) result of the multiplication.
"""
function inner(x::AbstractVector{N}, A::AbstractMatrix{N}, y::AbstractVector{N}
              ) where {N}
    dot(x, A * y)
end

"""
    _vector_type(T)

Return a corresponding vector type with respect to type `T`.

### Input

- `T` -- vector or matrix type

### Output

A vector type that corresponds in some sense (see Notes below) to `T`.

### Notes

- If `T` is a sparse vector or a sparse matrix, the corresponding type is also
  a sparse vector.
- If `T` is a regular vector (i.e. `Vector`) or a regular matrix (i.e. `Matrix`),
  the corresponding type is also a regular vector.
- Otherwise, the corresponding type is a regular vector.
"""
function _vector_type end

"""
    _matrix_type(T)

Return a corresponding matrix type with respect to type `T`.

### Input

- `T` -- vector type

### Output

A matrix type that corresponds in some sense (see Notes below) to `T`.

### Notes

- If `T` is a sparse vector or a sparse matrix, the corresponding type is also
  a sparse matrix.
- If `T` is a regular vector (i.e. `Vector`) or a regular matrix (i.e. `Matrix`),
  the corresponding type is also a regular matrix.
- Otherwise, the corresponding type is a regular matrix.
"""
function _matrix_type end

_vector_type(::Type{<:AbstractSparseArray{T}}) where T = SparseVector{T, Int}
_vector_type(VT::Type{<:AbstractVector{T}}) where T = VT
_vector_type(::Type{<:AbstractMatrix{T}}) where T = Vector{T}

_matrix_type(::Type{<:AbstractVector{T}}) where T = Matrix{T}
_matrix_type(MT::Type{<:AbstractMatrix{T}}) where T = MT
_matrix_type(::Type{<:AbstractSparseVector{T}}) where T = SparseMatrixCSC{T, Int}
