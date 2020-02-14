import Base: *, size, ∈, isempty

export SparseMatrixExp,
       ExponentialMap,
       size,
       ProjectionSparseMatrixExp,
       ExponentialProjectionMap,
       get_row,
       get_rows,
       get_column,
       get_columns

function load_expokit()
return quote

using .Expokit: expmv

end end  # quote / load_expokit


# --- SparseMatrixExp & ExponentialMap ---

"""
    SparseMatrixExp{N}

Type that represents the matrix exponential, ``\\exp(M)``, of a sparse matrix.

### Fields

- `M` -- sparse matrix; it should be square

### Examples

Take for example a random sparse matrix of dimensions ``100 × 100`` and with
occupation probability ``0.1``:

```jldoctest SparseMatrixExp_constructor
julia> using SparseArrays

julia> A = sprandn(100, 100, 0.1);

julia> using Expokit

julia> E = SparseMatrixExp(A);

julia> size(E)
(100, 100)
```

Here `E` is a lazy representation of ``\\exp(A)``. To compute with `E`, use
`get_row` and `get_column` resp. `get_rows` and `get_columns`. These functions
return row and column vectors (or matrices). For example:

```jldoctest SparseMatrixExp_constructor
julia> get_row(E, 10); # compute E[10, :]

julia> get_column(E, 10); # compute E[:, 10]

julia> get_rows(E, [10]); # same as get_row(E, 10) but a 1x100 matrix is returned

julia> get_columns(E, [10]); # same as get_column(E, 10) but a 100x1 matrix is returned
```

### Notes

This type is provided for use with very large and very sparse matrices.
The evaluation of the exponential matrix action over vectors relies on the
[Expokit](https://github.com/acroy/Expokit.jl) package. Hence, you will have to
install and load this optional dependency to have access to the functionality
of `SparseMatrixExp`.
"""
struct SparseMatrixExp{N, MN<:AbstractSparseMatrix{N}} <: AbstractMatrix{N}
    M::MN
    function SparseMatrixExp(M::MN) where {N, MN<:AbstractSparseMatrix{N}}
        @assert size(M, 1) == size(M, 2) "the lazy matrix exponential `SparseMatrixExp` " *
            "requires the given matrix to be square, but it has size $(size(M))"
        return new{N, MN}(M)
    end
end

SparseMatrixExp(M::Matrix) =
        error("only sparse matrices can be used to create a `SparseMatrixExp`")

Base.IndexStyle(::Type{<:SparseMatrixExp}) = IndexCartesian()
Base.getindex(spmexp::SparseMatrixExp, I::Vararg{Int, 2}) = get_column(spmexp, I[2])[I[1]]

function size(spmexp::SparseMatrixExp)
    return size(spmexp.M)
end

function size(spmexp::SparseMatrixExp, ax::Int)
    return size(spmexp.M, ax)
end

function get_column(spmexp::SparseMatrixExp{N}, j::Int) where {N}
    require(:Expokit; fun_name="get_column")

    n = size(spmexp, 1)
    aux = zeros(N, n)
    aux[j] = one(N)
    return expmv(one(N), spmexp.M, aux)
end

function get_columns(spmexp::SparseMatrixExp{N}, J::AbstractArray) where {N}
    require(:Expokit; fun_name="get_columns")

    n = size(spmexp, 1)
    aux = zeros(N, n)
    ans = zeros(N, n, length(J))
    count = 1
    one_N = one(N)
    zero_N = zero(N)
    @inbounds for j in J
        aux[j] = one_N
        ans[:, count] = expmv(one_N, spmexp.M, aux)
        aux[j] = zero_N
        count += 1
    end
    return ans
end

"""
    get_row(spmexp::SparseMatrixExp{N}, i::Int) where {N}

Return a single row of a sparse matrix exponential.

### Input

- `spmexp` -- sparse matrix exponential
- `i`      -- row index

### Output

A row vector corresponding to the `i`th row of the matrix exponential.

### Notes

This function uses Julia's `transpose` function to create the result.
The result is of type `Transpose`; in Julia versions older than v0.7, the result
was of type `RowVector`.
"""
function get_row(spmexp::SparseMatrixExp{N}, i::Int) where {N}
    require(:Expokit; fun_name="get_row")

    n = size(spmexp, 1)
    aux = zeros(N, n)
    aux[i] = one(N)
    return transpose(expmv(one(N), transpose(spmexp.M), aux))
end

function get_rows(spmexp::SparseMatrixExp{N}, I::AbstractArray{Int}) where {N}
    require(:Expokit; fun_name="get_rows")

    n = size(spmexp, 1)
    aux = zeros(N, n)
    ans = zeros(N, length(I), n)
    Mtranspose = transpose(spmexp.M)
    count = 1
    one_N = one(N)
    zero_N = zero(N)
    @inbounds for i in I
        aux[i] = one_N
        ans[count, :] = expmv(one_N, Mtranspose, aux)
        aux[i] = zero_N
        count += 1
    end
    return ans
end

"""
    ExponentialMap{N<:Real, S<:LazySet{N}} <: AbstractAffineMap{N, S}

Type that represents the action of an exponential map on a convex set.

### Fields

- `spmexp` -- sparse matrix exponential
- `X`      -- convex set

### Examples

The `ExponentialMap` type is overloaded to the usual times `*` operator when the
linear map is a lazy matrix exponential. For instance,

```jldoctest constructors
julia> using SparseArrays

julia> A = sprandn(100, 100, 0.1);

julia> E = SparseMatrixExp(A);

julia> B = BallInf(zeros(100), 1.);

julia> M = E * B; # represents the image set: exp(A) * B

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

julia> E * EmptySet(2)
EmptySet{Float64}(2)
```
"""
struct ExponentialMap{N<:Real, S<:LazySet{N}} <: AbstractAffineMap{N, S}
    spmexp::SparseMatrixExp{N}
    X::S
end

# ZeroSet is "almost absorbing" for ExponentialMap (only the dimension changes)
function ExponentialMap(spmexp::SparseMatrixExp{N}, Z::ZeroSet{N}
                       ) where {N<:Real}
    @assert dim(Z) == size(spmexp, 2) "an exponential map of size " *
            "$(size(spmexp)) cannot be applied to a set of dimension $(dim(Z))"
    return ZeroSet{N}(size(spmexp, 1))
end

isoperationtype(::Type{<:ExponentialMap}) = true
isconvextype(::Type{ExponentialMap{N, S}}) where {N, S} = isconvextype(S)

# EmptySet is absorbing for ExponentialMap
function ExponentialMap(spmexp::SparseMatrixExp{N}, ∅::EmptySet{N}
                       ) where {N<:Real}
    return ∅
end

"""
```
    *(spmexp::SparseMatrixExp{N}, X::LazySet{N}) where {N<:Real}
```

Return the exponential map of a convex set from a sparse matrix exponential.

### Input

- `spmexp` -- sparse matrix exponential
- `X`      -- convex set

### Output

The exponential map of the convex set.
"""
function *(spmexp::SparseMatrixExp{N}, X::LazySet{N}) where {N<:Real}
    return ExponentialMap(spmexp, X)
end


# --- AbstractAffineMap interface functions ---


function get_A(em::ExponentialMap)
    return em.spmexp
end

function get_b(em::ExponentialMap{N}) where {N<:Real}
    return spzeros(N, dim(em))
end

function get_X(em::ExponentialMap)
    return em.X
end


# --- LazySet interface functions ---


"""
    dim(em::ExponentialMap)

Return the dimension of an exponential map.

### Input

- `em` -- an ExponentialMap

### Output

The ambient dimension of the exponential map.
"""
function dim(em::ExponentialMap)
    return size(em.spmexp.M, 1)
end

"""
    σ(d::AbstractVector{N}, em::ExponentialMap{N}) where {N<:Real}

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
function σ(d::AbstractVector{N}, em::ExponentialMap{N}) where {N<:Real}
    require(:Expokit; fun_name="σ")

    d_dense = d isa Vector ? d : Vector(d)
    v = expmv(one(N), transpose(em.spmexp.M), d_dense) # v   <- exp(M') * d
    return expmv(one(N), em.spmexp.M, σ(v, em.X)) # res <- exp(M) * σ(v, S)
end

"""
    ρ(d::AbstractVector{N}, em::ExponentialMap{N}) where {N<:Real}

Return the support function of the exponential map.

### Input

- `d`  -- direction
- `em` -- exponential map

### Output

The support function in the given direction.

### Notes

If ``E = \\exp(M)⋅S``, where ``M`` is a matrix and ``S`` is a convex set, it
follows that ``ρ(d, E) = ρ(\\exp(M)^T d, S)`` for any direction ``d``.

We allow sparse direction vectors, but will convert them to dense vectors to be
able to use `expmv`.
"""
function ρ(d::AbstractVector{N}, em::ExponentialMap{N}) where {N<:Real}
    require(:Expokit; fun_name="ρ")

    d_dense = d isa Vector ? d : Vector(d)
    v = expmv(one(N), transpose(em.spmexp.M), d_dense) # v <- exp(M^T) * d
    return ρ(v, em.X)
end

"""
    ∈(x::AbstractVector{N}, em::ExponentialMap{N}) where {N<:Real}

Check whether a given point is contained in an exponential map of a convex set.

### Input

- `x`  -- point/vector
- `em` -- exponential map of a convex set

### Output

`true` iff ``x ∈ em``.

### Algorithm

This implementation exploits that ``x ∈ \\exp(M)⋅S`` iff ``\\exp(-M)⋅x ∈ S``.
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
function ∈(x::AbstractVector{N}, em::ExponentialMap{N}) where {N<:Real}
    require(:Expokit; fun_name="∈")

    @assert length(x) == dim(em)
    return expmv(-one(N), em.spmexp.M, x) ∈ em.X
end

"""
    vertices_list(em::ExponentialMap{N}) where {N<:Real}

Return the list of vertices of a (polytopic) exponential map.

### Input

- `em` -- exponential map

### Output

A list of vertices.

### Algorithm

We assume that the underlying set `X` is polytopic.
Then the result is just the exponential map applied to the vertices of `X`.
"""
function vertices_list(em::ExponentialMap{N}) where {N<:Real}
    require(:Expokit; fun_name="vertices_list")

    # collect low-dimensional vertices lists
    vlist_X = vertices_list(em.X)

    # create resulting vertices list
    vlist = Vector{Vector{N}}()
    sizehint!(vlist, length(vlist_X))
    for v in vlist_X
        push!(vlist, expmv(one(N), em.spmexp.M, v))
    end

    return vlist
end

"""
    isbounded(em::ExponentialMap)

Determine whether an exponential map is bounded.

### Input

- `em` -- exponential map

### Output

`true` iff the exponential map is bounded.
"""
function isbounded(em::ExponentialMap)
    return isbounded(em.X)
end

# --- ProjectionSparseMatrixExp & ExponentialProjectionMap ---

"""
    ProjectionSparseMatrixExp{N<:Real}

Type that represents the projection of a sparse matrix exponential, i.e.,
``L⋅\\exp(M)⋅R`` for a given sparse matrix ``M``.

### Fields

- `L` -- left multiplication matrix
- `E` -- sparse matrix exponential
- `R` -- right multiplication matrix
"""
struct ProjectionSparseMatrixExp{N<:Real, MN1<:AbstractSparseMatrix{N},
                                 MN2<:AbstractSparseMatrix{N},
                                 MN3<:AbstractSparseMatrix{N}}
    L::MN1
    spmexp::SparseMatrixExp{N, MN2}
    R::MN3
end

"""
    ExponentialProjectionMap{N<:Real, S<:LazySet{N}} <: LazySet{N}

Type that represents the application of a projection of a sparse matrix
exponential to a convex set.

### Fields

- `spmexp` -- projection of a sparse matrix exponential
- `X`      -- convex set
"""
struct ExponentialProjectionMap{N<:Real, S<:LazySet{N}} <: LazySet{N}
    projspmexp::ProjectionSparseMatrixExp
    X::S
end

isoperationtype(::Type{<:ExponentialProjectionMap}) = true
isconvextype(::Type{ExponentialProjectionMap{N, S}}) where {N, S} = isconvextype(S)

"""
```
    *(projspmexp::ProjectionSparseMatrixExp, X::LazySet)
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
function *(projspmexp::ProjectionSparseMatrixExp, X::LazySet)
    return ExponentialProjectionMap(projspmexp, X)
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
    σ(d::AbstractVector{N},
      eprojmap::ExponentialProjectionMap{N}) where {N<:Real}

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
           eprojmap::ExponentialProjectionMap{N}) where {N<:Real}
    require(:Expokit; fun_name="σ")

    d_dense = d isa Vector ? d : Vector(d)
    daux = transpose(eprojmap.projspmexp.L) * d_dense
    aux1 = expmv(one(N), transpose(eprojmap.projspmexp.spmexp.M), daux)
    daux = _At_mul_B(eprojmap.projspmexp.R, aux1)
    svec = σ(daux, eprojmap.X)

    aux2 = eprojmap.projspmexp.R * svec
    daux = expmv(one(N), eprojmap.projspmexp.spmexp.M, aux2)
    return eprojmap.projspmexp.L * daux
end

"""
    isbounded(eprojmap::ExponentialProjectionMap)

Determine whether an exponential projection map is bounded.

### Input

- `eprojmap` -- exponential projection map

### Output

`true` iff the exponential projection map is bounded.

### Algorithm

We first check if the left or right projection matrix is zero or the wrapped set
is bounded.
Otherwise, we check boundedness via [`isbounded_unit_dimensions`](@ref).
"""
function isbounded(eprojmap::ExponentialProjectionMap)
    if iszero(eprojmap.projspmexp.L) || iszero(eprojmap.projspmexp.R) ||
            isbounded(eprojmap.X)
        return true
    end
    return isbounded_unit_dimensions(eprojmap)
end

"""
    isempty(eprojmap::ExponentialProjectionMap)

Return if an exponential projection map is empty or not.

### Input

- `eprojmap` -- exponential projection map

### Output

`true` iff the wrapped set is empty.
"""
function isempty(eprojmap::ExponentialProjectionMap)
    return isempty(eprojmap.X)
end
