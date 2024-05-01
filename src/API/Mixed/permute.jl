"""
    permute(X::LazySet, p::AbstractVector{Int})

Permute the dimensions of a set according to a given permutation vector.

### Input

- `X` -- set
- `p` -- permutation vector

### Output

A new set corresponding to `X` where the dimensions have been permuted according
to `p`.
"""
function permute(::LazySet, ::AbstractVector{Int}) end
