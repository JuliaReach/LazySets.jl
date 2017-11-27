using Expokit

import Base: *, size

export SparseMatrixExp, ExponentialMap, size,
       ProjectionSparseMatrixExp, ExponentialProjectionMap,
       get_row, get_rows, get_column, get_columns

"""
    SparseMatrixExp

Type that represents the matrix exponential of a sparse matrix, and provides
evaluation of its action on vectors.

### Fields

- `M` -- sparse matrix

### Notes

This class is provided for use with very large and very sparse matrices. The
evaluation of the exponential matrix action over vectores relies on the
Expokit package.
"""
struct SparseMatrixExp{N<:Real}
    M::SparseMatrixCSC{N, Int64}

    # default constructor
    SparseMatrixExp{N}(M::SparseMatrixCSC{N, Int64}) where {N<:Real} = new{N}(M)
end
# type-less convenience constructor
SparseMatrixExp(M::SparseMatrixCSC{N, Int64}) where {N<:Real} =
    SparseMatrixExp{N}(M)

function size(spmexp::SparseMatrixExp)::Tuple{Int64,Int64}
    return size(spmexp.M)
end

function size(spmexp::SparseMatrixExp, ax::Int64)::Int64
    return size(spmexp.M, ax)
end

function get_column(spmexp::SparseMatrixExp, j::Int64)::Vector{Float64}
    n = size(spmexp, 1)
    aux = zeros(n)
    aux[j] = 1.0
    return expmv(1.0, spmexp.M, aux)
end

function get_columns(spmexp::SparseMatrixExp,
                     J::AbstractArray)::SparseMatrixCSC{Float64, Int64}
    n = size(spmexp, 1)
    aux = zeros(n)
    ans = spzeros(n, length(J))
    count = 1
    @inbounds for j in J
        aux[j] = 1.0
        ans[:, count] = sparse(expmv(1.0, spmexp.M, aux))
        aux[j] = 0.0
        count += 1
    end
    return ans
end

function get_row(spmexp::SparseMatrixExp, i::Int64)::Matrix{Float64}
    n = size(spmexp, 1)
    aux = zeros(n)
    aux[i] = 1.0
    return transpose(expmv(1.0, spmexp.M.', aux))
end

function get_rows(spmexp::SparseMatrixExp,
                  I::AbstractArray{Int64})::SparseMatrixCSC{Float64, Int64}
    n = size(spmexp, 1)
    aux = zeros(n)
    ans = spzeros(length(I), n)
    Mtranspose = spmexp.M.'
    count = 1
    @inbounds for i in I
        aux[i] = 1.0
        ans[count, :] = sparse(expmv(1.0, Mtranspose, aux))
        aux[i] = 0.0
        count += 1
    end
    return ans
end

"""
    ExponentialMap <: LazySet

Type that represents the action of an exponential map on a set.

### Fields

- `spmexp` -- a matrix exponential
- `X`      -- a convex set represented by its support function
"""
struct ExponentialMap{S<:LazySet} <: LazySet
    spmexp::SparseMatrixExp
    X::S
end
ExponentialMap(spmexp, X::S) where {S<:LazySet} = ExponentialMap{S}(spmexp,X)

# instantiate an exponential map from a sparse matrix exponential
function *(spmexp::SparseMatrixExp, X::LazySet)
    return ExponentialMap(spmexp, X)
end

"""
    dim(em)

The ambient dimension of a ExponentialMap.

### Input

- `em` -- an ExponentialMap
"""
function dim(em::ExponentialMap)::Int64
    return size(em.spmexp.M, 1)
end

# support vector of the exponential map
function σ(d::AbstractVector{Float64}, em::ExponentialMap)::AbstractVector{Float64}
    v = expmv(1.0, em.spmexp.M.', d)          # v   <- exp(A') * d
    res = expmv(1.0, em.spmexp.M, σ(v, em.X)) # res <- exp(A) * support_vector(v, S)
    return res
end

"""
    ProjectionSparseMatrixExp

Type that represents the projection of a SparseMatrixExp.

### Fields

- `L` -- left multiplication matrix
- `E` -- the exponential of a sparse matrix
- `R` -- right multiplication matrix

### Output

A type that abstract the matrix operation `L * exp(E.M) * R`, for a given sparse
matrix E.M.
"""
struct ProjectionSparseMatrixExp{N<:Real}
    L::SparseMatrixCSC{N, Int64}
    spmexp::SparseMatrixExp{N}
    R::SparseMatrixCSC{N, Int64}

    # default constructor
    ProjectionSparseMatrixExp{N}(L::SparseMatrixCSC{N, Int64},
                                 spmexp::SparseMatrixExp{N},
                                 R::SparseMatrixCSC{N, Int64}) where {N<:Real} =
        new{N}(L, spmexp, R)
end
# type-less convenience constructor
ProjectionSparseMatrixExp(L::SparseMatrixCSC{N, Int64},
                          spmexp::SparseMatrixExp{N},
                          R::SparseMatrixCSC{N, Int64}) where {N<:Real} =
    ProjectionSparseMatrixExp{N}(L, spmexp, R)

"""
    ExponentialProjectionMap

Type that represents the application of the projection of a SparseMatrixExp over
a given set.

### Fields

- `spmexp`   -- the projection of an exponential map
- `X`       -- a set represented by its support function
"""
struct ExponentialProjectionMap{S<:LazySet} <: LazySet
    projspmexp::ProjectionSparseMatrixExp
    X::S
end
ExponentialProjectionMap(projspmexp, X::S) where {S<:LazySet} =
    ExponentialProjectionMap{S}(projspmexp, X)

# instantiate an exponential map projection from matrix multiplication
function *(projspmexp::ProjectionSparseMatrixExp, X::LazySet)
    return ExponentialProjectionMap(projspmexp, X)
end

"""
    dim(eprojmap)

The ambient dimension of a ExponentialProjectionMap.

It is given by the output dimension (left-most matrix).

### Input

- `eprojmap` -- an ExponentialProjectionMap
"""
function dim(eprojmap::ExponentialProjectionMap)
    return size(eprojmap.projspmexp.L, 1)
end

"""
    σ(d, eprojmap)

Support vector of an `ExponentialProjectionMap`.

### Input

- `d`         -- a direction
- `eprojmap`  -- the projection of an exponential map

If `S = (LMR)B`, where `L` and `R` are dense matrices, `M` is a matrix
exponential, and `B` is a set, it follows that:
`σ(d, S) = LMR σ(R^T M^T L^T d, B)` for any direction `d`.
"""
function σ(d::AbstractVector{Float64},
           eprojmap::ExponentialProjectionMap)::Vector{Float64}
    daux = transpose(eprojmap.projspmexp.L) * d
    aux1 = expmv(1.0, eprojmap.projspmexp.spmexp.M.', daux)
    daux = transpose(eprojmap.projspmexp.R) * aux1
    aux1 = σ(daux, eprojmap.X)

    aux2 = eprojmap.projspmexp.R * aux1
    daux = expmv(1.0, eprojmap.projspmexp.spmexp.M, aux2)
    aux2 = eprojmap.projspmexp.L * daux
    return aux2
end
