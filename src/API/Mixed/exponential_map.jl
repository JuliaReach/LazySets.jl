"""
    exponential_map(M::AbstractMatrix, X::LazySet)

Compute the exponential map of `M` and `X`, i.e., ``eᴹ ⋅ X``.

### Input

- `M` -- matrix
- `X` -- set

### Output

A set representing the exponential map ``eᴹ ⋅ X``.
"""
function exponential_map(::AbstractMatrix, ::LazySet) end
