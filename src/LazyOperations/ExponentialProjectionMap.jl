export ProjectionSparseMatrixExp,
       ExponentialProjectionMap

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

function concretize(epm::ExponentialProjectionMap)
    p = epm.projspmexp
    M = p.L * exp(Matrix(p.spmexp.M)) * p.R
    return linear_map(M, concretize(epm.X))
end

"""
    dim(eprojmap::ExponentialProjectionMap)

Return the dimension of a projection of an exponential map.

### Input

- `eprojmap` -- projection of an exponential map

### Output

The ambient dimension of the projection of an exponential map.
"""
function dim(eprojmap::ExponentialProjectionMap)
    return size(eprojmap.projspmexp.L, 1)
end

"""
    σ(d::AbstractVector, eprojmap::ExponentialProjectionMap;
      [backend]=get_exponential_backend())

Return a support vector of a projection of an exponential map.

### Input

- `d`        -- direction
- `eprojmap` -- projection of an exponential map
- `backend`  -- (optional; default: `get_exponential_backend()`) exponentiation
                backend

### Output

A support vector in the given direction.
If the direction has norm zero, the result depends on the wrapped set.

### Notes

If ``S = (L⋅M⋅R)⋅X``, where ``L`` and ``R`` are matrices, ``M`` is a matrix
exponential, and ``X`` is a set, it follows that
``σ(d, S) = L⋅M⋅R⋅σ(R^T⋅M^T⋅L^T⋅d, X)`` for any direction ``d``.
"""
@validate function σ(d::AbstractVector, eprojmap::ExponentialProjectionMap;
                     backend=get_exponential_backend())
    N = promote_type(eltype(d), eltype(eprojmap))
    Lᵀd = transpose(eprojmap.projspmexp.L) * d
    eᴹLᵀd = _expmv(backend, one(N), transpose(eprojmap.projspmexp.spmexp.M), Lᵀd)
    RᵀeᴹLᵀd = At_mul_B(eprojmap.projspmexp.R, eᴹLᵀd)
    svec = σ(RᵀeᴹLᵀd, eprojmap.X)
    Rσ = eprojmap.projspmexp.R * svec
    eᴹRσ = _expmv(backend, one(N), eprojmap.projspmexp.spmexp.M, Rσ)
    return eprojmap.projspmexp.L * eᴹRσ
end

"""
    ρ(d::AbstractVector, eprojmap::ExponentialProjectionMap;
      [backend]=get_exponential_backend())

Evaluate the support function of a projection of an exponential map.

### Input

- `d`        -- direction
- `eprojmap` -- projection of an exponential map
- `backend`  -- (optional; default: `get_exponential_backend()`) exponentiation
                backend

### Output

Evaluation of the support function in the given direction.
If the direction has norm zero, the result depends on the wrapped set.

### Notes

If ``S = (L⋅M⋅R)⋅X``, where ``L`` and ``R`` are matrices, ``M`` is a matrix
exponential, and ``X`` is a set, it follows that ``ρ(d, S) = ρ(R^T⋅M^T⋅L^T⋅d, X)``
for any direction ``d``.
"""
@validate function ρ(d::AbstractVector, eprojmap::ExponentialProjectionMap;
                     backend=get_exponential_backend())
    N = promote_type(eltype(d), eltype(eprojmap))
    Lᵀd = transpose(eprojmap.projspmexp.L) * d
    eᴹLᵀd = _expmv(backend, one(N), transpose(eprojmap.projspmexp.spmexp.M), Lᵀd)
    RᵀeᴹLᵀd = At_mul_B(eprojmap.projspmexp.R, eᴹLᵀd)
    return ρ(RᵀeᴹLᵀd, eprojmap.X)
end

"""
    isbounded(eprojmap::ExponentialProjectionMap)

Check whether a projection of an exponential map is bounded.

### Input

- `eprojmap` -- projection of an exponential map

### Output

`true` iff the projection of an exponential map is bounded.

### Algorithm

We first check if the left or right projection matrix is zero or the wrapped set
is bounded.
Otherwise, we check boundedness via
[`LazySets._isbounded_unit_dimensions`](@ref).
"""
function isbounded(eprojmap::ExponentialProjectionMap)
    if iszero(eprojmap.projspmexp.L) || iszero(eprojmap.projspmexp.R) ||
       isbounded(eprojmap.X)
        return true
    end
    return _isbounded_unit_dimensions(eprojmap)
end
