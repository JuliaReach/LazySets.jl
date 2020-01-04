export _At_mul_B,
       _At_ldiv_B,
       DEFAULT_COND_TOL,
       hasfullrowrank,
       issquare,
       isinvertible,
       cross_product,
       delete_zero_columns!

# default tolerance for matrix condition number (see 'isinvertible')
const DEFAULT_COND_TOL = 1e6

# matrix-matrix multiplication
@inline _At_mul_B(A, B) = transpose(A) * B

# matrix-matrix division
@inline _At_ldiv_B(A, B) = transpose(A) \ B

@static if VERSION < v"1.2"
    # rank of sparse matrix (see JuliaLang #30415)
    LinearAlgebra.rank(M::SparseMatrixCSC) = rank(qr(M))
end

# rank of sparse submatrix (see #1497)
LinearAlgebra.rank(M::SubArray{N, 2, <:SparseMatrixCSC}) where {N} = rank(sparse(M))

"""
    issquare(M::AbstractMatrix)

Check whether a matrix is square.

### Input

- `M` -- matrix

### Output

`true` iff the matrix is square.
"""
function issquare(M::AbstractMatrix)
    m, n = size(M)
    return m == n
end

"""
    hasfullrowrank(M::AbstractMatrix)

Check whether a matrix has full row rank.

### Input

- `M` -- matrix

### Output

`true` iff the matrix has full row rank.
"""
function hasfullrowrank(M::AbstractMatrix)
    return rank(M) == size(M, 1)
end

"""
    isinvertible(M::Matrix; [cond_tol]::Number=DEFAULT_COND_TOL)

A sufficient check of a matrix being invertible (or nonsingular).

### Input

- `M`        -- matrix
- `cond_tol` -- (optional, default: `DEFAULT_COND_TOL`) tolerance of matrix
                condition

### Output

If the result is `true`, `M` is invertible.
If the result is `false`, the matrix is non-square or this function could not
conclude.

### Algorithm

We check whether the matrix is square and whether the
[matrix condition number](https://en.wikipedia.org/wiki/Condition_number#Matrices)
`cond(M)` is below some prescribed tolerance.
"""
function isinvertible(M::Matrix; cond_tol::Number=DEFAULT_COND_TOL)
    return issquare(M) && cond(M) < cond_tol
end

# cond is not available for sparse matrices; see JuliaLang#6485 and related issues
function isinvertible(M::AbstractSparseMatrix;
                      cond_tol::Number=DEFAULT_COND_TOL)
    return issquare(M) && isinvertible(Matrix(M), cond_tol=cond_tol)
end

function isinvertible(M::Diagonal; cond_tol=nothing)
    return !any(iszero, diag(M))
end

"""
    cross_product(M::AbstractMatrix{N}) where {N<:Real}

Compute the high-dimensional cross product of ``n-1`` ``n``-dimensional vectors.

### Input

- `M` -- ``n × n - 1``-dimensional matrix

### Output

A vector.

### Algorithm

The cross product is defined as follows:

```math
\\left[ \\dots, (-1)^{n+1} \\det(M^{[i]}), \\dots \\right]^T
```
where ``M^{[i]}`` is defined as ``M`` with the ``i``-th row removed.
See *Althoff, Stursberg, Buss: Computing Reachable Sets of Hybrid Systems Using
a Combination of Zonotopes and Polytopes. 2009.*
"""
function cross_product(M::AbstractMatrix{N}) where {N<:Real}
    n = size(M, 1)
    @assert 1 < n == size(M, 2) + 1 "the matrix must be n x (n-1) dimensional"

    v = Vector{N}(undef, n)
    minus = false
    for i in 1:n
        Mi = view(M, 1:n .!= i, :)  # remove i-th row
        d = det(Mi)
        if minus
            v[i] = -d
            minus = false
        else
            v[i] = d
            minus = true
        end
    end
    return v
end

# det cannot handle sparse matrices in some cases
cross_product(M::AbstractSparseMatrix) = cross_product(Matrix(M))
cross_product(M::SubArray{N, 2, <:AbstractSparseMatrix}) where {N} = cross_product(Matrix(M))

"""
    delete_zero_columns!(A::AbstractMatrix)

Remove all columns that only contain zeros from a given matrix.

### Input

- `A`    -- matrix
- `copy` -- (optional, default: `false`) flag to copy the matrix

### Output

A matrix.

If the input matrix `A` does not contain any zero column, we return `A` unless
the option `copy` is set.
If the input matrix contains zero columns, we always return a copy if the option
`copy` is set and otherwise a `SubArray` via `@view`.
"""
function delete_zero_columns!(A::AbstractMatrix, copy::Bool=false)
    n = size(A, 2)
    nonzero_columns = Vector{Int}()
    sizehint!(nonzero_columns, n)
    for i in 1:n
        if !iszero(A[:, i])
            push!(nonzero_columns, i)
        end
    end
    if copy
        if length(nonzero_columns) == n
            return copy(A)
        else
            return A[:, nonzero_columns]
        end
    else
        if length(nonzero_columns) == n
            return A
        else
            return @view A[:, nonzero_columns]
        end
    end
end
