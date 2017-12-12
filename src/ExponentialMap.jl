using Expokit

import Base: *, size

export SparseMatrixExp,
       ExponentialMap,
       size,
       ProjectionSparseMatrixExp,
       ExponentialProjectionMap,
       get_row,
       get_rows,
       get_column,
       get_columns

"""
    SparseMatrixExp{N<:Real}

Type that represents the matrix exponential, ``\\exp(M)``, of a sparse matrix.

### Fields

- `M` -- sparse matrix

### Notes

This type is provided for use with very large and very sparse matrices.
The evaluation of the exponential matrix action over vectors relies on the
[Expokit](https://github.com/acroy/Expokit.jl) package.
"""
struct SparseMatrixExp{N<:Real}
    M::SparseMatrixCSC{N, Int}
end

function size(spmexp::SparseMatrixExp)::Tuple{Int, Int}
    return size(spmexp.M)
end

function size(spmexp::SparseMatrixExp, ax::Int)::Int
    return size(spmexp.M, ax)
end

function get_column(spmexp::SparseMatrixExp, j::Int)::Vector{Float64}
    n = size(spmexp, 1)
    aux = zeros(n)
    aux[j] = 1.0
    return expmv(1.0, spmexp.M, aux)
end

function get_columns(spmexp::SparseMatrixExp,
                     J::AbstractArray)::SparseMatrixCSC{Float64, Int}
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

function get_row(spmexp::SparseMatrixExp, i::Int)::Matrix{Float64}
    n = size(spmexp, 1)
    aux = zeros(n)
    aux[i] = 1.0
    return transpose(expmv(1.0, spmexp.M.', aux))
end

function get_rows(spmexp::SparseMatrixExp,
                  I::AbstractArray{Int})::SparseMatrixCSC{Float64, Int}
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
    ExponentialMap{S<:LazySet} <: LazySet

Type that represents the action of an exponential map on a convex set.

### Fields

- `spmexp` -- sparse matrix exponential
- `X`      -- convex set
"""
struct ExponentialMap{S<:LazySet} <: LazySet
    spmexp::SparseMatrixExp
    X::S
end

"""
```
    *(spmexp::SparseMatrixExp, X::LazySet)::ExponentialMap
```

Return the exponential map of a convex set from a sparse matrix exponential.

### Input

- `spmexp` -- sparse matrix exponential
- `X`      -- convex set

### Output

The exponential map of the convex set.
"""
function *(spmexp::SparseMatrixExp, X::LazySet)::ExponentialMap
    return ExponentialMap(spmexp, X)
end

"""
    dim(em::ExponentialMap)::Int

Return the dimension of an exponential map.

### Input

- `em` -- an ExponentialMap

### Output

The ambient dimension of the exponential map.
"""
function dim(em::ExponentialMap)::Int
    return size(em.spmexp.M, 1)
end

"""
    σ(d::AbstractVector{Float64}, em::ExponentialMap)::AbstractVector{Float64}

Return the support vector of the exponential map.

### Input

- `d`  -- direction
- `em` -- exponential map

### Output

The support vector in the given direction.
If the direction has norm zero, the result depends on the wrapped set.

### Notes

If ``E = \\exp(M)⋅S``, where ``M`` is a matrix and ``S`` is a convex set, it
follows that ``σ(d, E) = \\exp(M)⋅σ(\\exp(M)^T d, S)`` for any direction ``d``.
"""
function σ(d::AbstractVector{Float64},
           em::ExponentialMap)::AbstractVector{Float64}
    v = expmv(1.0, em.spmexp.M.', d)           # v   <- exp(A') * d
    return expmv(1.0, em.spmexp.M, σ(v, em.X)) # res <- exp(A) * σ(v, S)
end

"""
    ProjectionSparseMatrixExp{N<:Real}

Type that represents the projection of a sparse matrix exponential, i.e.,
``L⋅\\exp(M)⋅R`` for a given sparse matrix ``M``.

### Fields

- `L` -- left multiplication matrix
- `E` -- sparse matrix exponential
- `R` -- right multiplication matrix
"""
struct ProjectionSparseMatrixExp{N<:Real}
    L::SparseMatrixCSC{N, Int}
    spmexp::SparseMatrixExp{N}
    R::SparseMatrixCSC{N, Int}
end

"""
    ExponentialProjectionMap{S<:LazySet} <: LazySet

Type that represents the application of a projection of a sparse matrix
exponential to a convex set.

### Fields

- `spmexp` -- projection of a sparse matrix exponential
- `X`      -- convex set
"""
struct ExponentialProjectionMap{S<:LazySet} <: LazySet
    projspmexp::ProjectionSparseMatrixExp
    X::S
end

"""
```
    *(projspmexp::ProjectionSparseMatrixExp, X::LazySet)::ExponentialProjectionMap
```

Return the application of a projection of a sparse matrix exponential to a
convex set.

### Input

- `projspmexp` -- projection of a sparse matrix exponential
- `X`          -- convex set

### Output

The application of the projection of a sparse matrix exponential to the convex
set.
"""
function *(projspmexp::ProjectionSparseMatrixExp,
           X::LazySet)::ExponentialProjectionMap
    return ExponentialProjectionMap(projspmexp, X)
end

"""
    dim(eprojmap::ExponentialProjectionMap)::Int

Return the dimension of a projection of an exponential map.

### Input

- `eprojmap` -- projection of an exponential map

### Output

The ambient dimension of the projection of an exponential map.
"""
function dim(eprojmap::ExponentialProjectionMap)::Int
    return size(eprojmap.projspmexp.L, 1)
end

"""
    σ(d::AbstractVector{Float64}, eprojmap::ExponentialProjectionMap)::AbstractVector{Float64}

Return the support vector of a projection of an exponential map.

### Input

- `d`        -- direction
- `eprojmap` -- projection of an exponential map

### Output

The support vector in the given direction.
If the direction has norm zero, the result depends on the wrapped set.

### Notes

If ``S = (L⋅M⋅R)⋅X``, where ``L`` and ``R`` are matrices, ``M`` is a matrix
exponential, and ``X`` is a set, it follows that
``σ(d, S) = L⋅M⋅R⋅σ(R^T⋅M^T⋅L^T⋅d, X)`` for any direction ``d``.
"""
function σ(d::AbstractVector{Float64},
           eprojmap::ExponentialProjectionMap)::AbstractVector{Float64}
    daux = transpose(eprojmap.projspmexp.L) * d
    aux1 = expmv(1.0, eprojmap.projspmexp.spmexp.M.', daux)
    daux = transpose(eprojmap.projspmexp.R) * aux1
    svec = σ(daux, eprojmap.X)

    aux2 = eprojmap.projspmexp.R * svec
    daux = expmv(1.0, eprojmap.projspmexp.spmexp.M, aux2)
    return eprojmap.projspmexp.L * daux
end
