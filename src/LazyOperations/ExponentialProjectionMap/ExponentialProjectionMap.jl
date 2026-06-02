"""
    ProjectionSparseMatrixExp{N, MN1<:AbstractSparseMatrix{N},
                                 MN2<:AbstractSparseMatrix{N},
                                 MN3<:AbstractSparseMatrix{N}}

Type that represents the projection of a sparse matrix exponential, i.e.,
``L⋅\\exp(M)⋅R`` for a given sparse matrix ``M``.

### Fields

- `L` -- left multiplication matrix
- `E` -- sparse matrix exponential
- `R` -- right multiplication matrix
"""
struct ProjectionSparseMatrixExp{N,MN1<:AbstractSparseMatrix{N},
                                 MN2<:AbstractSparseMatrix{N},
                                 MN3<:AbstractSparseMatrix{N}}
    L::MN1
    spmexp::SparseMatrixExp{N,MN2}
    R::MN3
end

"""
    ExponentialProjectionMap{N, S<:LazySet{N}} <: AbstractAffineMap{N, S}

Type that represents the application of a projection of a sparse matrix
exponential to a set.

### Fields

- `spmexp` -- projection of a sparse matrix exponential
- `X`      -- set

### Notes

The exponential projection preserves convexity: if `X` is convex, then any
exponential projection of `X` is convex as well.

### Examples

We use a random sparse projection matrix of dimensions ``6 × 6`` with occupation
probability ``0.5`` and apply it to the 2D unit ball in the infinity norm:

```jldoctest
julia> using SparseArrays

julia> R = sparse([5, 6], [1, 2], [1.0, 1.0]);

julia> L = sparse([1, 2], [1, 2], [1.0, 1.0], 2, 6);

julia> using ExponentialUtilities

julia> A = sprandn(6, 6, 0.5);

julia> E = SparseMatrixExp(A);

julia> M = ProjectionSparseMatrixExp(L, E, R);

julia> B = BallInf(zeros(2), 1.0);

julia> X = ExponentialProjectionMap(M, B);

julia> dim(X)
2
```
"""
struct ExponentialProjectionMap{N,S<:LazySet{N}} <: AbstractAffineMap{N,S}
    projspmexp::ProjectionSparseMatrixExp
    X::S
end

isoperationtype(::Type{<:ExponentialProjectionMap}) = true

isconvextype(::Type{ExponentialProjectionMap{N,S}}) where {N,S} = isconvextype(S)

"""
```
    *(projspmexp::ProjectionSparseMatrixExp, X::LazySet)
```

Alias to create an `ExponentialProjectionMap` object.

### Input

- `projspmexp` -- projection of a sparse matrix exponential
- `X`          -- set

### Output

The application of the projection of a sparse matrix exponential to the set.
"""
function *(projspmexp::ProjectionSparseMatrixExp, X::LazySet)
    return ExponentialProjectionMap(projspmexp, X)
end

function matrix(epm::ExponentialProjectionMap)
    projspmexp = epm.projspmexp
    return projspmexp.L * projspmexp.spmexp * projspmexp.R
end

function vector(epm::ExponentialProjectionMap{N}) where {N}
    return spzeros(N, dim(epm))
end

function set(epm::ExponentialProjectionMap)
    return epm.X
end

include("concretize.jl")
include("dim.jl")
include("isbounded.jl")
include("support_function.jl")
include("support_vector.jl")
