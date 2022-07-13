export QuadraticMap

"""
    QuadraticMap{N, S<:LazySet{N}, MVT<:AbstractVector{<:AbstractMatrix{N}}} <: LazySet{N}

Type that represents a quadratic map of a set.

### Fields

- `Q` -- matrices
- `X` -- set

### Notes

The quadratic map of set ``X`` for matrices ``Q\\^{(i)}`` is defined as

```math
\\left\\{ \\lambda | \\lambda_i = x^T Q\\^{(i)} x,~i = 1, \\ldots, n,~x \\in X \\right\\}
```

where each coordinate ``i`` is influenced by the ``i``-th matrix ``Q\\^{(i)}``.
"""
struct QuadraticMap{N, S<:LazySet{N}, MVT<:AbstractVector{<:AbstractMatrix{N}}}
    Q::MVT
    X::S

    # default constructor with dimension match check
    function QuadraticMap(Q::MVT, X::S) where {N, S<:LazySet{N}, MVT<:AbstractVector{<:AbstractMatrix{N}}}
        n = dim(X)
        @assert length(Q) == n "the number of matrices ($(length(Q))) needs " *
            "to match the dimension of the zonotope ($n)"
        @assert all(M -> checksquare(M) == n, Q) "dimension mismatch in the " *
            "matrices of the quadratic map applied to a set of size $n"
        return new{N, S, MVT}(Q, X)
    end
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
