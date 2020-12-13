export SingleEntryVector

"""
    SingleEntryVector{N} <: AbstractVector{N}

A sparse unit vector with arbitrary one-element.

### Fields

- `i` -- index of non-zero entry
- `n` -- vector length
- `v` -- non-zero entry
"""
struct SingleEntryVector{N} <: AbstractVector{N}
    i::Int
    n::Int
    v::N
end

# convenience constructor with one-element of corresponding type
SingleEntryVector{N}(i::Int, n::Int) where {N} =
    SingleEntryVector{N}(i, n, one(N))

function Base.getindex(e::SingleEntryVector{N}, i::Int) where {N}
    @boundscheck @assert 1 <= i <= e.n
    return i == e.i ? e.v : zero(N)
end

Base.size(e::SingleEntryVector) = (e.n,)

function Base.:(*)(A::AbstractMatrix{N}, e::SingleEntryVector{N}) where {N}
    return A[:, e.i] * e.v
end

function Base.:(*)(A::Transpose{N, <:AbstractMatrix{N}},
                   e::SingleEntryVector{N}) where {N}
    return A[:, e.i] * e.v
end

# multiplication with diagonal matrix
function Base.:(*)(D::Diagonal{N, V},
                   e::SingleEntryVector{N}) where {N, V<:AbstractVector{N}}
    return SingleEntryVector(e.i, e.n, D.diag[e.i] * e.v)
end

# negation
function Base.:(-)(e::SingleEntryVector{N}) where {N}
    return SingleEntryVector(e.i, e.n, -e.v)
end

function inner(e1::SingleEntryVector{N}, A::AbstractMatrix{N},
               e2::SingleEntryVector{N}) where {N}
    return A[e1.i, e2.i] * e1.v * e2.v
end

# norm
function LinearAlgebra.norm(e::SingleEntryVector, p::Real=Inf)
    return abs(e.v)
end
