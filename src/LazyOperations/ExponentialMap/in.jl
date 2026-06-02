"""
    in(x::AbstractVector, em::ExponentialMap;
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
@validate function in(x::AbstractVector, em::ExponentialMap;
                      backend=get_exponential_backend())
    N = promote_type(eltype(x), eltype(em))
    y = _expmv(backend, -one(N), em.expmat.M, x)
    return y ∈ em.X
end
