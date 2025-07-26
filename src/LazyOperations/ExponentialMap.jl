import Base: *, size

export SparseMatrixExp,
       ExponentialMap,
       size,
       get_row,
       get_rows,
       get_column,
       get_columns

"""
    SparseMatrixExp{N}

Type that represents the matrix exponential, ``\\exp(M)``, of a sparse matrix.

### Fields

- `M` -- sparse square matrix

### Examples

Take for example a random sparse matrix of dimensions ``100 × 100`` and with
occupation probability ``0.1``:

```jldoctest SparseMatrixExp_constructor
julia> using SparseArrays

julia> A = sprandn(100, 100, 0.1);

julia> using ExponentialUtilities

julia> E = SparseMatrixExp(A);

julia> size(E)
(100, 100)
```

Here `E` is a lazy representation of ``\\exp(A)``. To compute with `E`, use
`get_row` and `get_column` resp. `get_rows` and `get_columns`. These functions
return row and column vectors (or matrices). For example:

```jldoctest SparseMatrixExp_constructor
julia> get_row(E, 10);  # compute E[10, :]

julia> get_column(E, 10);  # compute E[:, 10]

julia> get_rows(E, [10]);  # same as get_row(E, 10), but yields a 1x100 matrix

julia> get_columns(E, [10]);  # same as get_column(E, 10), but yields a 100x1 matrix
```

### Notes

This type is provided for use with large and sparse matrices.
The evaluation of the exponential matrix action over vectors relies on external
packages such as
[ExponentialUtilities](https://github.com/SciML/ExponentialUtilities.jl) or
[Expokit](https://github.com/acroy/Expokit.jl).
Hence, you will have to install and load such an optional dependency to have
access to the functionality of `SparseMatrixExp`.
"""
struct SparseMatrixExp{N,MN<:AbstractSparseMatrix{N}} <: AbstractMatrix{N}
    M::MN

    # default constructor with dimension check
    function SparseMatrixExp(M::MN) where {N,MN<:AbstractSparseMatrix{N}}
        @assert size(M, 1) == size(M, 2) "the lazy matrix exponential " *
                                         "requires a square matrix, but it has size $(size(M))"
        return new{N,MN}(M)
    end
end

function SparseMatrixExp(M::AbstractMatrix)
    throw(ArgumentError("only sparse matrices can be used to create a `SparseMatrixExp`"))
end

Base.IndexStyle(::Type{<:SparseMatrixExp}) = IndexCartesian()

Base.getindex(spmexp::SparseMatrixExp, I::Vararg{Int,2}) = get_column(spmexp, I[2])[I[1]]

function size(spmexp::SparseMatrixExp)
    return size(spmexp.M)
end

function size(spmexp::SparseMatrixExp, ax::Int)
    return size(spmexp.M, ax)
end

"""
    get_column(spmexp::SparseMatrixExp{N}, j::Int;
               [backend]=get_exponential_backend()) where {N}

Compute a single column of a sparse matrix exponential.

### Input

- `spmexp`  -- sparse matrix exponential
- `j`       -- column index
- `backend` -- (optional; default: `get_exponential_backend()`) exponentiation
               backend

### Output

A column vector corresponding to the `j`th column of the matrix exponential.
"""
function get_column(spmexp::SparseMatrixExp{N}, j::Int;
                    backend=get_exponential_backend()) where {N}
    n = size(spmexp, 1)
    aux = zeros(N, n)
    aux[j] = one(N)
    return _expmv(backend, one(N), spmexp.M, aux)
end

"""
    get_columns(spmexp::SparseMatrixExp{N}, J::AbstractArray;
               [backend]=get_exponential_backend()) where {N}

Compute multiple columns of a sparse matrix exponential.

### Input

- `spmexp`  -- sparse matrix exponential
- `J`       -- list of column indices
- `backend` -- (optional; default: `get_exponential_backend()`) exponentiation
               backend

### Output

A matrix with one column vector for each entry in `J`.
"""
function get_columns(spmexp::SparseMatrixExp{N}, J::AbstractArray;
                     backend=get_exponential_backend()) where {N}
    n = size(spmexp, 1)
    aux = zeros(N, n)
    res = zeros(N, n, length(J))
    @inbounds for (k, j) in enumerate(J)
        aux[j] = one(N)
        res[:, k] = _expmv(backend, one(N), spmexp.M, aux)
        aux[j] = zero(N)
    end
    return res
end

"""
    get_row(spmexp::SparseMatrixExp{N}, i::Int;
            [backend]=get_exponential_backend()) where {N}

Compute a single row of a sparse matrix exponential.

### Input

- `spmexp`  -- sparse matrix exponential
- `i`       -- row index
- `backend` -- (optional; default: `get_exponential_backend()`) exponentiation
               backend

### Output

A row vector corresponding to the `i`th row of the matrix exponential.

### Notes

This implementation uses Julia's `transpose` function to create the result.
The result is of type `Transpose`; in Julia versions older than v0.7, the result
was of type `RowVector`.
"""
function get_row(spmexp::SparseMatrixExp{N}, i::Int;
                 backend=get_exponential_backend()) where {N}
    n = size(spmexp, 1)
    aux = zeros(N, n)
    aux[i] = one(N)
    return transpose(_expmv(backend, one(N), transpose(spmexp.M), aux))
end

"""
    get_rows(spmexp::SparseMatrixExp{N}, I::AbstractArray;
               [backend]=get_exponential_backend()) where {N}

Compute multiple rows of a sparse matrix exponential.

### Input

- `spmexp`  -- sparse matrix exponential
- `I`       -- list of row indices
- `backend` -- (optional; default: `get_exponential_backend()`) exponentiation
               backend

### Output

A transposed matrix with one (transposed) column vector for each entry in `I`.
"""
function get_rows(spmexp::SparseMatrixExp{N}, I::AbstractArray{Int};
                  backend=get_exponential_backend()) where {N}
    n = size(spmexp, 1)
    aux = zeros(N, n)
    res = zeros(N, length(I), n)
    Mtranspose = transpose(spmexp.M)
    @inbounds for (k, i) in enumerate(I)
        aux[i] = one(N)
        res[k, :] = _expmv(backend, one(N), Mtranspose, aux)
        aux[i] = zero(N)
    end
    return res
end

"""
    ExponentialMap{N,S<:LazySet{N},MAT<:Union{SparseMatrixExp{N},MatrixZonotopeExp{N}}} <:
       AbstractAffineMap{N,S}

Type that represents the action of an exponential map on a set.

### Fields

- `expmat` -- matrix exponential, either a concrete matrix exponential
              or a matrix zonotope exponential represented by an `MatrixZonotopeExp`
- `X`      -- set

### Notes

The exponential map preserves convexity: if `X` is convex, then any exponential
map of `X` is convex as well.

### Examples

The `ExponentialMap` type is overloaded to the usual times (`*`) operator when
the linear map is a lazy matrix exponential. For instance:

```jldoctest constructors
julia> using SparseArrays

julia> A = sprandn(100, 100, 0.1);

julia> E = SparseMatrixExp(A);

julia> B = BallInf(zeros(100), 1.);

julia> M = E * B;  # represents the set: exp(A) * B

julia> M isa ExponentialMap
true

julia> dim(M)
100
```

The application of an `ExponentialMap` to a `ZeroSet` or an `EmptySet` is
simplified automatically.

```jldoctest constructors
julia> E * ZeroSet(100)
ZeroSet{Float64}(100)

julia> E * EmptySet(100)
∅(100)
```
"""
struct ExponentialMap{N,S<:LazySet{N},MAT<:Union{SparseMatrixExp{N},MatrixZonotopeExp{N}}} <:
       AbstractAffineMap{N,S}
    expmat::MAT
    X::S
    function ExponentialMap(expmat::MAT,
                            X::S) where {N,S<:LazySet{N},NM,
                                         MAT<:Union{SparseMatrixExp{NM},MatrixZonotopeExp{NM}}}
        @assert dim(X) == size(expmat, 2) "an exponential map of size $(size(expmat)) cannot " *
                                          "be applied to a set of dimension $(dim(X))"
        @assert size(expmat, 1) == size(expmat, 2) "the matrix exponential requires a square " *
                                                   "matrix, but it has size $(size(expmat))"
        return new{N,S,MAT}(expmat, X)
    end
end

# ZeroSet is "almost absorbing" for ExponentialMap (only the dimension changes)
function ExponentialMap(spmexp::SparseMatrixExp, Z::ZeroSet)
    N = promote_type(eltype(spmexp), eltype(Z))
    @assert dim(Z) == size(spmexp, 2) "an exponential map of size " *
                                      "$(size(spmexp)) cannot be applied to a set of dimension $(dim(Z))"
    return ZeroSet{N}(size(spmexp, 1))
end

# EmptySet is "almost absorbing" for ExponentialMap (only the dimension changes)
function ExponentialMap(spmexp::SparseMatrixExp, ∅::EmptySet)
    N = promote_type(eltype(spmexp), eltype(∅))
    @assert dim(∅) == size(spmexp, 2) "an exponential map of size " *
                                      "$(size(spmexp)) cannot be applied to a set of dimension $(dim(∅))"
    return EmptySet{N}(size(spmexp, 1))
end

# auto-convert M to SparseMatrixExp
function ExponentialMap(M::AbstractMatrix, X::LazySet)
    return ExponentialMap(SparseMatrixExp(M), X)
end

isoperationtype(::Type{<:ExponentialMap}) = true

isconvextype(::Type{ExponentialMap{N,S}}) where {N,S} = isconvextype(S)

"""
```
    *(spmexp::SparseMatrixExp, X::LazySet)
```

Alias to create an `ExponentialMap` object.

### Input

- `spmexp` -- sparse matrix exponential
- `X`      -- set

### Output

The exponential map of the set.
"""
function *(spmexp::SparseMatrixExp, X::LazySet)
    return ExponentialMap(spmexp, X)
end

function matrix(em::ExponentialMap)
    return em.expmat
end

function vector(em::ExponentialMap{N}) where {N}
    return spzeros(N, dim(em))
end

function set(em::ExponentialMap)
    return em.X
end

"""
    dim(em::ExponentialMap)

Return the dimension of an exponential map.

### Input

- `em` -- exponential map

### Output

The ambient dimension of the exponential map.
"""
function dim(em::ExponentialMap)
    return size(em.expmat, 1)
end

"""
    σ(d::AbstractVector, em::ExponentialMap;
      [backend]=get_exponential_backend())

Return a support vector of an exponential map.

### Input

- `d`       -- direction
- `em`      -- exponential map
- `backend` -- (optional; default: `get_exponential_backend()`) exponentiation
               backend

### Output

A support vector in the given direction.
If the direction has norm zero, the result depends on the wrapped set.

### Notes

If ``E = \\exp(M)⋅X``, where ``M`` is a matrix and ``X`` is a set, it
follows that ``σ(d, E) = \\exp(M)⋅σ(\\exp(M)^T d, X)`` for any direction ``d``.
"""
@validate function σ(d::AbstractVector, em::ExponentialMap;
                     backend=get_exponential_backend())
    N = promote_type(eltype(d), eltype(em))
    v = _expmv(backend, one(N), transpose(em.expmat.M), d)  # exp(M^T) * d
    return _expmv(backend, one(N), em.expmat.M, σ(v, em.X))  # exp(M) * σ(v, X)
end

"""
    ρ(d::AbstractVector, em::ExponentialMap;
      [backend]=get_exponential_backend())

Evaluate the support function of the exponential map.

### Input

- `d`       -- direction
- `em`      -- exponential map
- `backend` -- (optional; default: `get_exponential_backend()`) exponentiation
               backend

### Output

The evaluation of the support function in the given direction.

### Notes

If ``E = \\exp(M)⋅X``, where ``M`` is a matrix and ``X`` is a set, it
follows that ``ρ(d, E) = ρ(\\exp(M)^T d, X)`` for any direction ``d``.
"""
@validate function ρ(d::AbstractVector, em::ExponentialMap;
                     backend=get_exponential_backend())
    N = promote_type(eltype(d), eltype(em))
    v = _expmv(backend, one(N), transpose(em.expmat.M), d)  # exp(M^T) * d
    return ρ(v, em.X)
end

function concretize(em::ExponentialMap)
    return exponential_map(Matrix(em.expmat.M), concretize(em.X))
end

"""
    ∈(x::AbstractVector, em::ExponentialMap;
      [backend]=get_exponential_backend())

Check whether a given point is contained in an exponential map of a set.

### Input

- `x`       -- point/vector
- `em`      -- exponential map of a set
- `backend` -- (optional; default: `get_exponential_backend()`) exponentiation
               backend

### Output

`true` iff ``x ∈ em``.

### Algorithm

This implementation exploits that ``x ∈ \\exp(M)⋅X`` iff ``\\exp(-M)⋅x ∈ X``.
This follows from ``\\exp(-M)⋅\\exp(M) = I`` for any ``M``.

### Examples

```jldoctest
julia> using SparseArrays

julia> em = ExponentialMap(
        SparseMatrixExp(sparse([1, 2], [1, 2], [2.0, 1.0], 2, 2)),
        BallInf([1., 1.], 1.));

julia> [-1.0, 1.0] ∈ em
false
julia> [1.0, 1.0] ∈ em
true
```
"""
@validate function ∈(x::AbstractVector, em::ExponentialMap;
                     backend=get_exponential_backend())
    N = promote_type(eltype(x), eltype(em))
    y = _expmv(backend, -one(N), em.expmat.M, x)
    return y ∈ em.X
end

"""
    vertices_list(em::ExponentialMap; [backend]=get_exponential_backend())

Return the list of vertices of a (polytopic) exponential map.

### Input

- `em`      -- polytopic exponential map
- `backend` -- (optional; default: `get_exponential_backend()`) exponentiation
               backend

### Output

A list of vertices.

### Algorithm

We assume that the underlying set `X` is polytopic.
Then the result is just the exponential map applied to the vertices of `X`.
"""
function vertices_list(em::ExponentialMap; backend=get_exponential_backend())
    # collect vertices lists of wrapped set
    vlist_X = vertices_list(em.X)

    # create resulting vertices list
    N = eltype(em)
    vlist = Vector{Vector{N}}(undef, length(vlist_X))
    @inbounds for (i, v) in enumerate(vlist_X)
        vlist[i] = _expmv(backend, one(N), em.expmat.M, v)
    end

    return vlist
end

"""
    isbounded(em::ExponentialMap)

Check whether an exponential map is bounded.

### Input

- `em` -- exponential map

### Output

`true` iff the exponential map is bounded.
"""
function isbounded(em::ExponentialMap)
    return isbounded(em.X)
end

function isboundedtype(::Type{<:ExponentialMap{N,S}}) where {N,S}
    return isboundedtype(S)
end
