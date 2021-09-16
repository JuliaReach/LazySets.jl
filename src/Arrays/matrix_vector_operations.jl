export inner,
       _abs_sum,
       to_matrix

# computes ‖a^T G‖₁
@inline function _abs_sum(a::AbstractVector, G::AbstractMatrix)
   n, p = size(G)
   N = promote_type(eltype(a), eltype(G))
   abs_sum = zero(N)
   @inbounds for j in 1:p
       aux = zero(N)
       @simd for i in 1:n
           aux += a[i] * G[i, j]
       end
       abs_sum += abs(aux)
   end
   return abs_sum
end

# computes ‖a^T G‖₁ for `a` being a sparse vector
@inline function _abs_sum(a::AbstractSparseVector, G::AbstractMatrix)
   return sum(abs, transpose(a) * G)
end

# computes ‖a^T G‖₁ for `a` having only one nonzero element
@inline function _abs_sum(a::SingleEntryVector, G::AbstractMatrix)
   p = size(G, 2)
   i = a.i
   v = abs(a.v)
   N = promote_type(eltype(a), eltype(G))
   abs_sum = zero(N)
   @inbounds for j in 1:p
       abs_sum += abs(G[i, j])
   end
   abs_sum *= v
   return abs_sum
end

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

- `T` -- vector or matrix type

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

# matrix constructors
_matrix(m, n, MT::Type{<:AbstractMatrix{T}}) where T = Matrix{T}(undef, m, n)
_matrix(m, n, MT::Type{<:SparseMatrixCSC{T}}) where T = spzeros(T, m, n)

"""
    to_matrix(vectors::AbstractVector{VN},
              [m]=length(vectors[1]),
              [MT]=_matrix_type(VN)) where {VN}

### Input

- `vectors` -- list of vectors
- `m`       -- (optional; default: `length(vectors[1])`) number of rows
- `MT`      -- (optional; default: `_matrix_type(VN)`) type of target matrix

### Output

A matrix with the column vectors from `vectors` in the same order.
"""
function to_matrix(vectors::AbstractVector{VN},
                   m=length(vectors[1]),
                   mat_type=_matrix_type(VN)) where {VN}
    n = length(vectors)
    M = _matrix(m, n, mat_type)
    @inbounds for (j, vj) in enumerate(vectors)
        M[:, j] = vj
    end
    return M
end

# no-op
_similar_type(x::AbstractVector) = typeof(x)

function load_copy_finalize_static()

return quote
    _similar_type(x::StaticArrays.StaticArray) = StaticArrays.similar_type(x)
end # quote

end # end load_copy_finalize_static
