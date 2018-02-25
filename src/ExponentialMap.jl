using Expokit

import Base: *, size, ∈

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
    SparseMatrixExp{N}

Type that represents the matrix exponential, ``\\exp(M)``, of a sparse matrix.

### Fields

- `M` -- sparse matrix

### Notes

This type is provided for use with very large and very sparse matrices.
The evaluation of the exponential matrix action over vectors relies on the
[Expokit](https://github.com/acroy/Expokit.jl) package.
"""
struct SparseMatrixExp{N}
    M::SparseMatrixCSC{N, Int}
end

function size(spmexp::SparseMatrixExp)::Tuple{Int, Int}
    return size(spmexp.M)
end

function size(spmexp::SparseMatrixExp, ax::Int)::Int
    return size(spmexp.M, ax)
end

function get_column(spmexp::SparseMatrixExp{N}, j::Int)::Vector{N} where {N}
    n = size(spmexp, 1)
    aux = zeros(N, n)
    aux[j] = one(N)
    return expmv(one(N), spmexp.M, aux)
end

function get_columns(spmexp::SparseMatrixExp{N},
                     J::AbstractArray)::SparseMatrixCSC{N, Int} where {N}
    n = size(spmexp, 1)
    aux = zeros(N, n)
    ans = spzeros(N, n, length(J))
    count = 1
    one_N = one(N)
    zero_N = zero(N)
    @inbounds for j in J
        aux[j] = one_N
        ans[:, count] = sparse(expmv(one_N, spmexp.M, aux))
        aux[j] = zero_N
        count += 1
    end
    return ans
end

function get_row(spmexp::SparseMatrixExp{N}, i::Int)::Matrix{N} where {N}
    n = size(spmexp, 1)
    aux = zeros(N, n)
    aux[i] = one(N)
    return transpose(expmv(one(N), spmexp.M.', aux))
end

function get_rows(spmexp::SparseMatrixExp{N},
                  I::AbstractArray{Int})::SparseMatrixCSC{N, Int} where {N}
    n = size(spmexp, 1)
    aux = zeros(N, n)
    ans = spzeros(N, length(I), n)
    Mtranspose = spmexp.M.'
    count = 1
    one_N = one(N)
    zero_N = zero(N)
    @inbounds for i in I
        aux[i] = one_N
        ans[count, :] = sparse(expmv(one_N, Mtranspose, aux))
        aux[i] = zero_N
        count += 1
    end
    return ans
end

"""
    ExponentialMap{N<:Real, S<:LazySet{N}} <: LazySet{N}

Type that represents the action of an exponential map on a set.

### Fields

- `spmexp` -- sparse matrix exponential
- `X`      -- set
"""
struct ExponentialMap{N, S<:LazySet{N}} <: LazySet{N}
    spmexp::SparseMatrixExp{N}
    X::S
end
# type-less convenience constructor
ExponentialMap(spmexp::SparseMatrixExp, X::S) where {S<:LazySet{N}} where {N} =
    ExponentialMap{N, S}(spmexp, X)

"""
```
    *(spmexp::SparseMatrixExp, X::LazySet)::ExponentialMap
```

Return the exponential map of a set from a sparse matrix exponential.

### Input

- `spmexp` -- sparse matrix exponential
- `X`      -- set

### Output

The exponential map of the set.
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
    σ(d::AbstractVector{N}, em::ExponentialMap)::AbstractVector{N} where {N<:Real}

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

We allow sparse direction vectors, but will convert them to dense vectors to be
able to use `expmv`.
"""
function σ(d::AbstractVector{N}, em::ExponentialMap)::AbstractVector{N} where {N<:Real}
    d_dense = d isa Vector ? d : Vector(d)
    v = expmv(one(N), em.spmexp.M.', d_dense)     # v   <- exp(A') * d
    return expmv(one(N), em.spmexp.M, σ(v, em.X)) # res <- exp(A) * σ(v, S)
end

"""
    ∈(x::AbstractVector{N}, em::ExponentialMap{<:LazySet{N}})::Bool where {N<:Real}

Check whether a given point is contained in an exponential map of a set.

### Input

- `x`  -- point/vector
- `em` -- exponential map of a set

### Output

`true` iff ``x ∈ em``.

### Algorithm

This implementation exploits that ``x ∈ \\exp(M)⋅S`` iff ``\\exp(-M)⋅x ∈ S``.
This follows from ``\\exp(-M)⋅\\exp(M) = I`` for any ``M``.

### Examples

```jldoctest
julia> em = ExponentialMap(SparseMatrixExp(SparseMatrixCSC([2.0 0.0; 0.0 1.0])),
                           BallInf([1., 1.], 1.));

julia> ∈([5.0, 1.0], em)
false
julia> ∈([1.0, 1.0], em)
true
```
"""
function ∈(x::AbstractVector{N}, em::ExponentialMap{N, <:LazySet{N}})::Bool where {N<:Real}
    @assert length(x) == dim(em)
    return ∈(exp.(-em.spmexp.M) * x, em.X)
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
    ExponentialProjectionMap{N<:Real, S<:LazySet{N}} <: LazySet{N}

Type that represents the application of a projection of a sparse matrix
exponential to a set.

### Fields

- `spmexp` -- projection of a sparse matrix exponential
- `X`      -- set
"""
struct ExponentialProjectionMap{N<:Real, S<:LazySet{N}} <: LazySet{N}
    projspmexp::ProjectionSparseMatrixExp
    X::S
end
# type-less convenience constructor
ExponentialProjectionMap(projspmexp::ProjectionSparseMatrixExp, X::S
                        ) where {S<:LazySet{N}} where {N<:Real} =
    ExponentialProjectionMap{N, S}(projspmexp, X)

"""
```
    *(projspmexp::ProjectionSparseMatrixExp,
      X::LazySet)::ExponentialProjectionMap
```

Return the application of a projection of a sparse matrix exponential to a set.

### Input

- `projspmexp` -- projection of a sparse matrix exponential
- `X`          -- set

### Output

The application of the projection of a sparse matrix exponential to the set.
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
    σ(d::AbstractVector{N},
      eprojmap::ExponentialProjectionMap)::AbstractVector{N} where {N<:Real}

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

We allow sparse direction vectors, but will convert them to dense vectors to be
able to use `expmv`.
"""
function σ(d::AbstractVector{N},
           eprojmap::ExponentialProjectionMap)::AbstractVector{N} where {N<:Real}
    d_dense = d isa Vector ? d : Vector(d)
    daux = transpose(eprojmap.projspmexp.L) * d_dense
    aux1 = expmv(one(N), eprojmap.projspmexp.spmexp.M.', daux)
    daux = transpose(eprojmap.projspmexp.R) * aux1
    svec = σ(daux, eprojmap.X)

    aux2 = eprojmap.projspmexp.R * svec
    daux = expmv(one(N), eprojmap.projspmexp.spmexp.M, aux2)
    return eprojmap.projspmexp.L * daux
end
