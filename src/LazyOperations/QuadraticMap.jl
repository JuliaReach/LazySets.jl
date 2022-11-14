export QuadraticMap

"""
    QuadraticMap{N, S<:LazySet{N}, MVT<:AbstractVector{<:AbstractMatrix{N}}}
        <: LazySet{N}

Type that represents a quadratic map of a set.

### Fields

- `Q` -- matrices
- `X` -- set

### Notes

The quadratic map of a set ``X`` given ``n`` square matrices ``Q^{(i)}`` is
defined as

```math
\\left\\{ \\lambda \\mid \\lambda_i = x^T Q^{(i)} x,~i = 1, \\ldots, n,~x \\in X \\right\\}
```

where each coordinate ``i`` is influenced by the ``i``-th matrix ``Q^{(i)}``.
"""
struct QuadraticMap{N, S<:LazySet{N}, MVT<:AbstractVector{<:AbstractMatrix{N}}} <: LazySet{N}
    Q::MVT
    X::S

    # default constructor with dimension match check
    function QuadraticMap(Q::MVT, X::S) where {N, S<:LazySet{N},
                                       MVT<:AbstractVector{<:AbstractMatrix{N}}}
        n = dim(X)
        @assert length(Q) == n "the number of matrices ($(length(Q))) needs " *
            "to match the dimension of the set ($n)"
        @assert all(M -> checksquare(M) == n, Q) "dimension mismatch in the " *
            "matrices of the quadratic map applied to a set of dimension $n"
        return new{N, S, MVT}(Q, X)
    end
end

function isoperationtype(P::Type{<:QuadraticMap})
    return true
end

function isconvextype(P::Type{<:QuadraticMap})
    return false
end

function isboundedtype(P::Type{<:QuadraticMap{MVT, S}}) where {MVT, S}
    return isboundedtype(S)
end

"""
    dim(qm::QuadraticMap)

Return the dimension of a quadratic map.

### Input

- `qm` -- quadratic map

### Output

The ambient dimension of the quadratic map.
"""
function dim(qm::QuadraticMap)
    return dim(qm.X)
end
